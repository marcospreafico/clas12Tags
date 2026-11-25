#!/usr/bin/env python3
"""Extract the eight active GenericTrap definitions from the pinned BDX source."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'experiments/bdx/geometry_crs60.pl'


def arithmetic(value):
    def calculate(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return float(node.value)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Div, ast.Mult, ast.Add, ast.Sub)):
            a, b = calculate(node.left), calculate(node.right)
            if isinstance(node.op, ast.Div): return a / b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Add): return a + b
            return a - b
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub): return -calculate(node.operand)
        raise ValueError(f'Unsupported source expression: {value}')
    return calculate(ast.parse(value.strip(), mode='eval').body)


def load_types(source=SOURCE):
    text = source.read_text()
    active = re.sub(r'#[^\n]*', '', text)
    match = re.search(r'my\s+%bgo_types\s*=\s*\((.*?)\n\);', active, re.S)
    if not match or 'G4GenericTrap' not in active:
        raise ValueError('Expected BDX GenericTrap definitions not found')
    types = {}
    for key, body in re.findall(r'"([1-8])"\s*=>\s*\{(.*?)\}', match[1], re.S):
        if key in types: raise ValueError(f'Duplicate crystal {key}')
        values = {k: arithmetic(v) for k, v in re.findall(r'(\w+)\s*=>\s*([^,]+),', body)}
        if set(values) != {'Bg', 'Bm', 'H', 'bg', 'bm', 'h', 'L', 'theta'}:
            raise ValueError(f'Unexpected dimensions for crystal {key}')
        # Source quantities are half dimensions in cm. theta is unused by its
        # active G4GenericTrap construction: do not introduce an extra rotation.
        dims = {k: v * 10 for k, v in values.items() if k != 'theta'}
        bm, bg, h, Bm, Bg, H = (dims[k] for k in ('bm','bg','h','Bm','Bg','H'))
        vertices = [[-bm,-h],[-bg,h],[bg,h],[bm,-h],[-Bm,-H],[-Bg,H],[Bg,H],[Bm,-H]]
        types[int(key)] = dict(half_dimensions_mm=dims, vertices_mm=vertices,
            half_length_mm=dims['L'], source_theta_deg_unused=values['theta'],
            # Area centroid, not the centre of the bounding rectangle.
            sensor_center_mm=[0, H*(Bg-Bm)/(3*(Bg+Bm)), dims['L']])
    if set(types) != set(range(1,9)): raise ValueError('Expected exactly eight active crystal types')
    return types


def geometry_text(shape):
    d = shape['half_dimensions_mm']
    dimensions = ' '.join(f'{x:.12g}*mm' for x in [d['L'], *sum(shape['vertices_mm'], [])])
    # Extrude the exact large-base quadrilateral for 0.1 mm beyond z=+L.
    face = shape['vertices_mm'][4:]
    collector = ' '.join(f'{x:.12g}*mm' for x in [0.05, *sum(face + face, [])])
    def row(name, description, pos, solid, dims, material, color):
        return ' | '.join([name,'root',description,pos,'0*deg 0*deg 0*deg',color,solid,
            dims,material,'no','1','1','1','1','1','no','no','no']) + '\n'
    return row('bdx_bgo','Isolated BDX BGO crystal','0*mm 0*mm 0*mm',
        'G4GenericTrap',dimensions,'G4_BGO','55bbff') + row('bdx_flux',
        'First-arrival optical flux on entire large base',f'0*mm 0*mm {d["L"]+0.05:.12g}*mm',
        'G4GenericTrap',collector,'G4_Galactic','ff3333')


def card_text(crystal_type):
    card = ET.Element('gcard')
    ET.SubElement(card,'detector', name='bdx_optics', factory='TEXT', variation=f'bgo{crystal_type}')
    for key, value in dict(PHYSICS='FTFP_BERT + STD + Optical', HALL_MATERIAL='G4_AIR',
        HALL_DIMENSIONS='300*mm,300*mm,300*mm', BEAM_P='mu-,4*GeV,90*deg,0*deg',
        BEAM_V='(-100,0,0)mm', SPREAD_P='0*GeV,0*deg,0*deg', SPREAD_V='(0,0)mm',
        USE_GUI='0', OUTPUT='hipo,optical.hipo', FILTER_HITS='0', N='1', PRINT_EVENT='1', ENERGY_CUT='0',
        PRODUCTIONCUT='0.1*mm', RANDOM='12345').items():
        ET.SubElement(card,'option',name=key,value=value)
    ET.indent(card)
    return ET.tostring(card,encoding='unicode') + '\n'


def generate(output, source=SOURCE):
    output.mkdir(parents=True,exist_ok=True)
    types = load_types(source)
    manifest = dict(source=str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else str(source),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        upstream_commit='1d2313f8cc701e43bf461bee44ea12f43c7d402d',
        units='mm; lengths are half dimensions unless stated otherwise',
        collector='full_large_base_first_arrival', schema_version=2, types=types)
    for key, shape in types.items():
        (output/f'bdx_optics__geometry_bgo{key}.txt').write_text(geometry_text(shape))
        (output/f'bgo{key}.gcard').write_text(card_text(key))
        for table in ('materials','hit','parameters','mirrors'):
            (output/f'bdx_optics__{table}_bgo{key}.txt').write_text('')
    (output/'bdx_optics__bank.txt').write_text('')
    (output/'geometry_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'experiments/bdx_optics')
    args=parser.parse_args()
    manifest=generate(args.output)
    print(f'Generated eight crystal variations in {args.output}; source SHA256 {manifest["source_sha256"]}')
