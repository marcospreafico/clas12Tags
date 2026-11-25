#!/usr/bin/env python3
"""One small independent optical job; publish only after HIPO validation."""
import json,math,os,subprocess,sys,tarfile
from pathlib import Path
from release import sha,write_json


def execute(command,log,env=None,cwd=None):
    with log.open('w') as stream:
        result=subprocess.run(command,cwd=cwd,env=env,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:
        print(log.read_text(errors='replace')[-12000:],file=sys.stderr)
        raise RuntimeError('Job stage failed: '+str(log))


def run(jobfile):
    jobfile=Path(jobfile).absolute(); job=json.loads(jobfile.read_text())
    runtime=Path(os.environ['BDX_RUNTIME_ROOT']).resolve()
    work=Path.cwd()/'run';work.mkdir(exist_ok=False)
    for name,text in job['files'].items():
        if Path(name).name!=name: raise ValueError('Invalid input name')
        (work/name).write_text(text)
    (work/'optics.conf').write_text(''.join(f'{k}={v}\n' for k,v in job['configuration'].items()))
    env=dict(os.environ,BDX_OPTICS_CONFIG=str(work/'optics.conf'),BDX_OPTICS_CSV=str(work/'events.csv'),
             BDX_OPTICS_TYPE=str(job['crystal_type']),BDX_OPTICS_SEED=str(job['seed']))
    env.pop('BDX_OPTICS_ARRIVALS',None)
    command=[str(runtime/'source/gemc'),'run.gcard',f'-N={job["events"]}',f'-RANDOM={job["seed"]}',
             f'-BEAM_P={job["particle"]},{job["momentum_GeV"]}*GeV,90*deg,0*deg',
             f'-BEAM_V=(-100,0,{120-job["distance_mm"]:.12g})mm','-USE_GUI=0','-OUTPUT=hipo,optical.hipo']
    execute(command,work/'gemc.log',env,work)
    execute([str(runtime/'build/bdx_optics/bdx_hipo'),str(work/'optical.hipo'),'--reduce',str(work/'reduction.json')],work/'validation.log')
    reduced=json.loads((work/'reduction.json').read_text())['events']
    if len(reduced)!=job['events']: raise ValueError('HIPO event count mismatch')
    for row in reduced:
        c=row['config']; e=row['event']
        for name,value in [('crystal_type',job['crystal_type']),('seed',job['seed']),('collector',1)]:
            if c[name]!=value: raise ValueError('HIPO metadata mismatch: '+name)
        for name,value in job['configuration'].items():
            expected=({'paint':1,'contact':2,'bare':0}[value] if name=='wrap_model' else float(value))
            if not math.isclose(c[name],expected,rel_tol=1e-10,abs_tol=1e-10): raise ValueError('HIPO optical setting mismatch: '+name)
        if not math.isclose(e['nominal_distance_mm'],job['distance_mm'],abs_tol=1e-8): raise ValueError('HIPO source distance mismatch')
    receipt=dict(format='BDX_SWIF_RESULT_V1',status='completed',job=job['name'],crystal_type=job['crystal_type'],
        distance_mm=job['distance_mm'],events=len(reduced),seed=job['seed'],particle=job['particle'],
        momentum_GeV=job['momentum_GeV'],runtime_sha256=job['runtime_sha256'],
        job_sha256=sha(jobfile),hipo_sha256=sha(work/'optical.hipo'),
        generated_scint=sum(r['event']['generated_scint'] for r in reduced),
        arrived_scint=sum(r['selected_scint'] for r in reduced))
    write_json(work/'result.json',receipt)
    with tarfile.open(work/'logs.tar.gz','w:gz') as archive:
        for name in ['gemc.log','validation.log','events.csv','optics.conf','run.gcard']:
            archive.add(work/name,arcname=name)
    (work/'reduction.json').unlink()
    print(json.dumps(receipt))

if __name__=='__main__': run(sys.argv[1])
