#!/usr/bin/env python3
"""Aggregate many BDX HIPO jobs, retaining event denominators and checking provenance."""
import argparse,csv,glob,hashlib,json,math,os,subprocess,tempfile
from pathlib import Path
from bdx_geometry import ROOT
from bdx_data import IMAGE
from analyze_bdx_optics import summarize,svg_plot


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()


def canonical(value): return json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':'))


def discover(args):
    files=[]; expected={}; missing=[]
    for campaign in args.campaign:
        campaign=campaign.resolve(); plan=json.loads((campaign/'plan.json').read_text())
        for job in plan['jobs']:
            path=campaign/job['hipo']; receipt=campaign/job['receipt']
            if not path.is_file() or not receipt.is_file():
                missing.append(str(path)); continue
            expected[path.resolve()]=dict(job=job,runtime_sha256=plan['runtime_sha256'])
            files.append(path.resolve())
    if missing and not args.allow_partial: raise ValueError(f'{len(missing)} planned results missing; wait or explicitly use --allow-partial')
    patterns=list(args.inputs)
    for listing in args.input_list:
        for line in listing.read_text().splitlines():
            line=line.strip()
            if line and not line.startswith('#'): patterns.append(str(listing.resolve().parent/line))
    for pattern in patterns:
        found=sorted(glob.glob(pattern,recursive=True))
        if not found: raise ValueError('No files matched: '+pattern)
        files.extend(Path(p).resolve() for p in found)
    if not files: raise ValueError('No HIPO files')
    if len(set(files))!=len(files): raise ValueError('Duplicate input paths (possibly overlapping patterns/campaigns)')
    return files,expected,missing


def reduce_file(path,destination,args):
    reader=(args.reader or ROOT/'build/bdx_optics/bdx_hipo').resolve()
    if not reader.is_file(): raise ValueError('Compile the updated bdx_hipo reader or supply --reader')
    options=['--reduce',str(destination)]
    if args.square_mm is not None: options.append(str(args.square_mm))
    if args.native: command=[str(reader),str(path),*options]
    else:
        # One file at a time; large photon banks are reduced in C++, never exported to CSV.
        options=['--reduce','/bdx-output/'+destination.name]
        if args.square_mm is not None: options.append(str(args.square_mm))
        command=['docker','run','--rm','--platform','linux/amd64','--entrypoint','/bin/bash',
            '-v',f'{reader}:/bdx-reader:ro','-v',f'{path.parent}:/bdx-input:ro',
            '-v',f'{destination.parent}:/bdx-output',os.environ.get('BDX_GEMC_IMAGE',IMAGE),'-lc',
            'module load hipo/4.5.0; exec "$@"','bdx-reduce','/bdx-reader','/bdx-input/'+path.name,*options]
    with destination.with_suffix('.log').open('w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
    result=json.loads(destination.read_text())
    if result['format']!='BDX_OPTICS_REDUCTION_V1' or not result['events']: raise ValueError('Missing/unsupported optical events')
    return result['events']


class Aggregator:
    def __init__(self):
        self.groups={}; self.geometry={}; self.signature=None; self.file_seeds=set(); self.sources=[]; self.hashes=set(); self.runtime=None

    def add(self,path,digest,rows,receipt=None,expected=None):
        if digest in self.hashes: raise ValueError('Duplicate HIPO content: '+str(path))
        self.hashes.add(digest)
        if receipt:
            if receipt['status']!='completed' or receipt['hipo_sha256']!=digest: raise ValueError('Invalid receipt/hash: '+str(path))
            if len(rows)!=receipt['events']: raise ValueError('Receipt event count mismatch')
            if self.runtime is not None and self.runtime!=receipt['runtime_sha256']: raise ValueError('Mixing native build releases is not supported')
            self.runtime=receipt['runtime_sha256']
        if expected:
            if not receipt or receipt['runtime_sha256']!=expected['runtime_sha256']: raise ValueError('Result differs from campaign build')
            job=expected['job']
            if any(receipt[k]!=job[k] for k in ['seed','events','crystal_type','distance_mm']): raise ValueError('Result differs from planned job')
        identity=None; point=None
        for row in rows:
            cfg=row['config']; event=dict(row['event']); primaries=row['primaries']
            kind=int(cfg['crystal_type']); distance=float(event['nominal_distance_mm'])
            if cfg['schema_version']!=2 or cfg['collector']!=1 or not 1<=kind<=8: raise ValueError('Requires schema-2 full-base BGO production')
            if receipt and (cfg['seed']!=receipt['seed'] or kind!=receipt['crystal_type'] or not math.isclose(distance,receipt['distance_mm'],abs_tol=1e-8)): raise ValueError('HIPO config differs from receipt')
            dimensions={k:cfg[k] for k in ['half_length_mm','Bg_mm','Bm_mm','H_mm','bg_mm','bm_mm','h_mm','center_y_mm','plane_z_mm']}
            if kind in self.geometry and self.geometry[kind]!=dimensions: raise ValueError('Mixed geometries for one crystal type')
            self.geometry[kind]=dimensions
            optics={k:v for k,v in cfg.items() if k not in dimensions and k not in ['crystal_type','seed']}
            if not primaries: raise ValueError('Missing primary particles')
            beam=[{k:v for k,v in particle.items() if k!='vz'} for particle in primaries]
            signature=canonical(dict(optics=optics,beam=beam))
            if self.signature is not None and signature!=self.signature: raise ValueError('Incompatible optical settings or primary beam; analyze these datasets separately')
            self.signature=signature
            current=(kind,distance,canonical(cfg),canonical(beam))
            if identity is not None and current!=identity: raise ValueError('One input file must contain a single type, distance and seed')
            identity=current; point=(kind,distance)
            generated=sum(int(event[k]) for k in ['generated_scint','generated_cherenkov','injected_optical'])
            terminated=sum(int(event[k]) for k in ['arrived_scint','arrived_cherenkov','arrived_injected','bulk_absorbed','surface_absorbed','escaped'])
            if generated!=terminated: raise ValueError('Optical budget mismatch')
            selected=int(row['selected_scint'])
            if not 0<=selected<=event['arrived_scint']: raise ValueError('Invalid selected count')
            event['arrived_scint']=selected
            # Only event summaries are retained; photon rows never enter Python.
            self.groups.setdefault(point,[]).append(event)
        if identity in self.file_seeds: raise ValueError('Repeated seed/configuration at one scan point: '+str(path))
        self.file_seeds.add(identity)
        self.sources.append(dict(path=str(path),sha256=digest,events=len(rows),crystal_type=point[0],distance_mm=point[1],
            receipt=bool(receipt),seed=rows[0]['config']['seed']))

    def results(self):
        config=json.loads(self.signature)['optics']
        data=[dict(crystal_type=k[0],distance_mm=k[1],files=sum(s['crystal_type']==k[0] and s['distance_mm']==k[1] for s in self.sources),
                   **summarize(rows,float(config['light_yield_per_MeV']))) for k,rows in sorted(self.groups.items())]
        return data,config


def analyze(args):
    if args.square_mm is not None and (not math.isfinite(args.square_mm) or args.square_mm<=0): raise ValueError('Invalid square size')
    files,expected,missing=discover(args)
    output=args.output.resolve()
    # Protect raw data and campaign inputs from accidental analysis writes.
    for campaign in args.campaign:
        if output.is_relative_to(campaign.resolve()): raise ValueError('Analysis output must be outside production campaign')
    output.mkdir(parents=True,exist_ok=False); accumulator=Aggregator()
    for index,path in enumerate(files):
        digest=sha(path); receipt_path=path.with_suffix('.json')
        receipt=json.loads(receipt_path.read_text()) if receipt_path.is_file() else None
        if receipt and receipt.get('format')!='BDX_SWIF_RESULT_V1': receipt=None
        with tempfile.TemporaryDirectory(prefix='reduce_',dir=output) as directory:
            rows=reduce_file(path,Path(directory)/'events.json',args)
        accumulator.add(path,digest,rows,receipt,expected.get(path))
        print(f'[{index+1}/{len(files)}] {path.name}: {len(rows)} events',flush=True)
    data,config=accumulator.results()
    with (output/'collection.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    selection=dict(area='full_base' if args.square_mm is None else 'square',square_side_mm=args.square_mm,origin='scintillation')
    # Per-point event counts can differ, e.g. incomplete campaigns or final short jobs.
    manifest=dict(configuration=config,events_per_point='variabili',geometry=dict(collector='full_large_base_first_arrival'))
    plotting_version=svg_plot(data,output/'collection.svg',manifest,selection)
    reader=(args.reader or ROOT/'build/bdx_optics/bdx_hipo').resolve()
    report=dict(format='BDX_MULTI_HIPO_ANALYSIS_V1',partial=bool(missing),missing=missing,
        selection=selection,configuration=config,geometry=accumulator.geometry,
        runtime_sha256=accumulator.runtime,sources=accumulator.sources,points=data,
        normalisation='ratio of total selected scintillation arrivals to total generated scintillation photons',
        uncertainty='event-cluster ratio SEM; Wilson68 is photon counting precision only',
        analysis_sha256=sha(Path(__file__)),statistics_sha256=sha(ROOT/'scripts/analyze_bdx_optics.py'),
        reader_sha256=sha(reader),plotting_version=plotting_version)
    (output/'analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (output/'inputs.list').write_text(''.join(str(path)+'\n' for path in files))
    print(f'Aggregated {len(files)} files into {len(data)} scan points: {output}')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',nargs='*',help='HIPO files or quoted glob patterns')
    p.add_argument('--input-list',type=Path,action='append',default=[])
    p.add_argument('--campaign',type=Path,action='append',default=[])
    p.add_argument('--allow-partial',action='store_true',help='explicitly include only completed planned jobs; mark analysis partial')
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--square-mm',type=float)
    p.add_argument('--native',action='store_true');p.add_argument('--reader',type=Path)
    a=p.parse_args(); analyze(a)

if __name__=='__main__':
    try: main()
    except (OSError,ValueError,RuntimeError,subprocess.SubprocessError) as error: raise SystemExit(str(error))
