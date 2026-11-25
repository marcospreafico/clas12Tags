#!/usr/bin/env python3
"""Run reproducible GEMC scans of isolated BDX BGO crystals (Docker or native)."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from bdx_geometry import ROOT, SOURCE, generate
from bdx_data import IMAGE, hipo_command




def digest(path):
    checksum=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024),b''): checksum.update(block)
    return checksum.hexdigest()


def config_values(path, overrides):
    values = {}
    for line in path.read_text().splitlines():
        line = line.split('#',1)[0].strip()
        if not line: continue
        key, value = map(str.strip,line.split('=',1))
        if key in values: raise ValueError(f'Duplicate config key {key}')
        values[key] = value
    allowed={'light_yield_per_MeV','refractive_index','decay_time_ns','wavelength_nm',
             'absorption_length_mm','wrap_reflectivity','wrap_model','cherenkov',
             'diffuse_length_mm','transition_length_mm','specular_length_mm','boundary_push_mm'}
    if set(values) != allowed: raise ValueError(f'Config keys must be exactly {sorted(allowed)}')
    for override in overrides:
        key, value = override.split('=',1)
        if key not in allowed: raise ValueError(f'Unknown setting {key}')
        values[key] = value
    for key in allowed-{'wrap_model'}:
        if not math.isfinite(float(values[key])): raise ValueError(f'Nonfinite {key}')
    if float(values['light_yield_per_MeV']) <= 0: raise ValueError('Light yield must be positive')
    if not 0<=float(values['boundary_push_mm'])<=1e-3: raise ValueError('Invalid numerical boundary push')
    if float(values['refractive_index']) < 1: raise ValueError('Index must be >=1')
    if float(values['decay_time_ns']) <= 0: raise ValueError('Decay time must be positive')
    if not 250 < float(values['wavelength_nm']) < 800: raise ValueError('Wavelength outside table support')
    if float(values['absorption_length_mm']) < 0: raise ValueError('Absorption length must be nonnegative')
    if not 0 <= float(values['wrap_reflectivity']) <= 1: raise ValueError('Reflectivity outside [0,1]')
    if values['wrap_model'] not in {'paint','contact','bare'}: raise ValueError('Unknown wrapping model')
    lengths=[float(values[k]) for k in ('diffuse_length_mm','transition_length_mm','specular_length_mm')]
    if any(x<0 for x in lengths) or abs(sum(lengths)-240)>1e-6:
        raise ValueError('Nonnegative coating lengths must sum to 240 mm')
    if values['cherenkov'] not in {'0','1'}: raise ValueError('cherenkov must be 0 or 1')
    if (values['wrap_model'] != 'bare' and float(values['wrap_reflectivity']) == 1
            and float(values['absorption_length_mm']) == 0):
        raise ValueError('Perfect lossless reflector can trap rays indefinitely')
    return values


def atomic_json(path, data):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(path)


def validate_events(path, requested):
    with path.open() as stream: rows=list(csv.DictReader(stream))
    if len(rows) != requested: raise RuntimeError(f'Expected {requested} events, found {len(rows)}')
    ids=[]
    for row in rows:
        ids.append(int(row['event']))
        generated=sum(int(row[k]) for k in ('generated_scint','generated_cherenkov','injected_optical'))
        terminated=sum(int(row[k]) for k in ('arrived_scint','arrived_cherenkov','arrived_injected',
                                             'bulk_absorbed','surface_absorbed','escaped'))
        if generated != terminated: raise RuntimeError('Photon budget failed')
    if ids != list(range(requested)): raise RuntimeError('Unexpected event sequence')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--types',nargs='+',type=int,default=list(range(1,9)),choices=range(1,9))
    parser.add_argument('--distances-mm',nargs='+',type=float,default=[5,25,50,75,100,125,150,175,200,225,235])
    parser.add_argument('--events',type=int,default=10,help='events per type and distance')
    parser.add_argument('--seed',type=int,default=20261123)
    parser.add_argument('--momentum-GeV',type=float,default=4)
    parser.add_argument('--particle',default='mu-',choices=['mu-','mu+','e-','e+','proton','pi-','pi+','gamma'])
    parser.add_argument('--config',type=Path,default=ROOT/'experiments/bdx_optics/geometry_only.conf')
    parser.add_argument('--set',action='append',default=[],metavar='KEY=VALUE',help='override optical settings')
    parser.add_argument('--output',type=Path,help='new directory inside the clone (never overwritten)')
    parser.add_argument('--native',action='store_true',help='use local source/gemc with loaded GEMC environment')
    parser.add_argument('--dry-run',action='store_true',help='prepare all inputs without running GEMC')
    args=parser.parse_args(argv)
    if args.events <= 0 or not 0 < args.seed < 2147483647 or not math.isfinite(args.momentum_GeV) or args.momentum_GeV <= 0:
        parser.error('Positive events, seed in [1,2147483646] and positive finite momentum required')
    if len(set(args.types)) != len(args.types): parser.error('Duplicate crystal types')
    if len(set(args.distances_mm)) != len(args.distances_mm): parser.error('Duplicate distances')
    if any(not math.isfinite(d) or not 0 < d < 240 for d in args.distances_mm):
        parser.error('Distances must be inside the 240 mm crystal: 0 < d < 240')
    values=config_values(args.config,args.set)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output=(args.output or ROOT/'build/bdx_optics'/stamp).resolve()
    if not output.is_relative_to(ROOT): parser.error('Output must be inside the clone for the Docker mount')
    executable=ROOT/'source/gemc'
    if not args.dry_run and not executable.is_file():
        parser.error('First compile GEMC using scripts/build_bdx_optics.sh')
    if not args.dry_run and not (ROOT/'build/bdx_optics/bdx_hipo').is_file():
        parser.error('Compile GEMC and its HIPO reader using scripts/build_bdx_optics.sh')
    output.mkdir(parents=True,exist_ok=False)
    generated=output/'inputs'
    geometry=generate(generated)
    shutil.copy2(SOURCE,generated/'geometry_crs60.source.pl')
    config=generated/'optics.conf'
    config.write_text(''.join(f'{k}={v}\n' for k,v in values.items()))
    image=os.environ.get('BDX_GEMC_IMAGE',IMAGE)
    commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
    files=[SOURCE,ROOT/'source/src/BdxCrystalOptics.cc',ROOT/'source/src/BdxCrystalOptics.h',
           ROOT/'source/src/MDetectorConstruction.cc',ROOT/'source/src/MSteppingAction.cc',
           ROOT/'source/src/MEventAction.cc',ROOT/'source/physics/PhysicsList.cc',
           Path(__file__),ROOT/'scripts/bdx_geometry.py',ROOT/'scripts/bdx_data.py',ROOT/'scripts/bdx_hipo.cc',ROOT/'source/src/BdxOpticalSchema.h',
           ROOT/'source/output/hipoSchemas.cc',ROOT/'source/output/hipo_output.h',ROOT/'source/output/hipo_output.cc',
           ROOT/'scripts/build_bdx_optics.sh',ROOT/'source/SConstruct',
           ROOT/'source/output/outputFactory.h',ROOT/'source/output/outputFactory.cc']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in files}
    manifest=dict(format='BDX_OPTICS_CAMPAIGN_V2', created_utc=stamp, upstream_commit=commit,
        geometry=geometry, configuration=values, source_sha256=hashes,
        gemc_sha256=digest(executable) if executable.is_file() else None,
        git_status=subprocess.check_output(['git','-C',str(ROOT),'status','--short'],text=True),
        image=image, particle=args.particle, momentum_GeV=args.momentum_GeV,
        seed=args.seed, events_per_point=args.events, status='planned', runs=[])
    if not args.native and not args.dry_run:
        manifest['image_id']=subprocess.check_output(['docker','image','inspect',image,'--format','{{.Id}}'],text=True).strip()
    plan=output/'campaign.json'
    seeds=set()
    for kind in args.types:
        for distance in args.distances_mm:
            key=f'{args.seed}:{kind}:{distance:.12g}'
            seed=int.from_bytes(hashlib.sha256(key.encode()).digest()[:8],'big') % 2147483646 + 1
            if seed in seeds: raise ValueError('Derived seed collision; choose another base seed')
            seeds.add(seed)
            run=output/f'bgo{kind}'/f'd{distance:.12g}mm'
            run.mkdir(parents=True)
            for p in generated.glob('bdx_optics__*'):
                # Snapshot inputs required by this variation, plus the bank table.
                if f'_bgo{kind}.' in p.name or p.name.endswith('__bank.txt'): shutil.copy2(p,run/p.name)
            shutil.copy2(generated/f'bgo{kind}.gcard',run/'run.gcard')
            shutil.copy2(config,run/'optics.conf')
            z=geometry['types'][kind]['half_length_mm']-distance
            options=['run.gcard',f'-N={args.events}',f'-RANDOM={seed}',
                     f'-BEAM_P={args.particle},{args.momentum_GeV}*GeV,90*deg,0*deg',
                     f'-BEAM_V=(-100,0,{z:.12g})mm']
            if args.native:
                command=[str(executable),*options]
            else:
                inside='/work/clas12Tags/'+str(run.relative_to(ROOT))
                command=['docker','run','--rm','--platform','linux/amd64','--entrypoint','/bin/bash',
                    '-v',f'{ROOT}:/work/clas12Tags','-w',inside,
                    '-e','BDX_OPTICS_CONFIG='+inside+'/optics.conf',
                    '-e','BDX_OPTICS_CSV='+inside+'/events.csv',
                    '-e',f'BDX_OPTICS_TYPE={kind}','-e',f'BDX_OPTICS_SEED={seed}']
                command += [image,'-lc',
                    'module load ccdb/1.99.6 hipo/4.5.0; '
                    'export GEMC=/work/clas12Tags GEMC_DATA_DIR=/work/clas12Tags; '
                    'exec /work/clas12Tags/source/gemc "$@"','bdx-run',*options]
            record=dict(crystal_type=kind,distance_mm=distance,seed=seed,
                        directory=str(run.relative_to(output)),status='planned',command=command)
            manifest['runs'].append(record)
    atomic_json(plan,manifest)
    if args.dry_run:
        print(f'Prepared {len(manifest["runs"])} points: {plan}')
        return 0
    manifest['status']='running'
    for record in manifest['runs']:
        run=output/record['directory']
        record['status']='running'
        atomic_json(plan,manifest)
        print(f'BGO {record["crystal_type"]}, d={record["distance_mm"]} mm, seed={record["seed"]}',flush=True)
        environment=dict(os.environ)
        environment.pop('BDX_OPTICS_ARRIVALS',None)
        if args.native:
            environment.update(BDX_OPTICS_CONFIG=str(run/'optics.conf'),BDX_OPTICS_CSV=str(run/'events.csv'))
            environment.update(BDX_OPTICS_TYPE=str(record['crystal_type']),BDX_OPTICS_SEED=str(record['seed']))
        failure_log=run/'gemc.log'
        try:
            with failure_log.open('w') as log:
                subprocess.run(record['command'],cwd=run,env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
            validate_events(run/'events.csv',args.events)
            failure_log=run/'hipo_check.log'
            with failure_log.open('w') as log:
                subprocess.run(hipo_command(run/'optical.hipo','--check',run/'banks.json',args.native,image),
                               check=True,stdout=log,stderr=subprocess.STDOUT)
            inventory=json.loads((run/'banks.json').read_text())
            if inventory['events']!=args.events: raise RuntimeError('Missing HIPO events')
            record['hipo_file']='optical.hipo'
            record['hipo_sha256']=digest(run/'optical.hipo')
            record['bank_rows']=inventory['banks']
            record['events_sha256']=digest(run/'events.csv')
            record['status']='completed'
        except (subprocess.CalledProcessError,RuntimeError,OSError) as error:
            record['status']='failed';record['error']=str(error);manifest['status']='failed'
            atomic_json(plan,manifest)
            record['failure_log']=failure_log.name
            atomic_json(plan,manifest)
            print(f'Failed: {failure_log}',file=sys.stderr)
            try:
                lines=failure_log.read_text(errors='replace').splitlines()
                reason=next((line.strip() for line in reversed(lines) if 'what():' in line),None)
                if reason:
                    record['failure_reason']=reason
                    atomic_json(plan,manifest)
                    print(reason,file=sys.stderr)
            except OSError:
                pass
            return 1
        atomic_json(plan,manifest)
    manifest['status']='completed'
    atomic_json(plan,manifest)
    print(f'Production completed (no analysis): {output}')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError,FileExistsError) as error:
        raise SystemExit(str(error))
