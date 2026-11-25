#!/usr/bin/env python3
"""Prepare or explicitly submit a frozen multi-job JLab SWIF2 campaign."""
import argparse,getpass,hashlib,json,math,re,shlex,shutil,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from bdx_geometry import generate,SOURCE
from run_bdx_optics import config_values
from release import sha,write_json


def prepare(a):
    if not a.account or not re.fullmatch(r'[\w.-]+',a.account): raise ValueError('Supply your assigned SWIF2 --account')
    if not a.release or not a.output: raise ValueError('Supply --release and a new shared --output directory')
    release_path=a.release.resolve(); build=json.loads(release_path.read_text())
    archive=Path(build['archive'])
    if sha(archive)!=build['archive_sha256']: raise ValueError('Build archive changed')
    # Allow dirty development checkouts: the complete source was captured at build.
    for name,digest in build['source']['source_sha256'].items():
        if name=='farm/ifarm.env': continue
        if name.endswith(('.py','.sh','.cc','.h','.pl')) or name.endswith('SConstruct'):
            if not (ROOT/name).is_file() or sha(ROOT/name)!=digest: raise ValueError('Source differs from native build; rebuild: '+name)
    if a.events_per_job<1 or a.jobs_per_point<1 or (a.events_per_point is not None and a.events_per_point<1): raise ValueError('Event/job counts must be positive')
    if not 0<a.seed<2147483647: raise ValueError('Seed must be in [1,2147483646]')
    if len(a.types)!=len(set(a.types)) or len(a.distances_mm)!=len(set(a.distances_mm)): raise ValueError('Repeated type or distance')
    if any(not math.isfinite(d) or not 0<d<240 for d in a.distances_mm): raise ValueError('Distances must be inside (0,240) mm')
    if not math.isfinite(a.momentum_GeV) or a.momentum_GeV<=0: raise ValueError('Invalid momentum')
    output=a.output.resolve()
    if not re.fullmatch(r'[/A-Za-z0-9_.-]+',str(output)): raise ValueError('Use a shared output path without whitespace/shell metacharacters')
    workflow=a.workflow or 'bdx_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',workflow): raise ValueError('Invalid workflow name')
    cfg=config_values(a.config,a.set)
    output.mkdir(parents=True,exist_ok=False); inputs=output/'inputs'; inputs.mkdir(); (output/'results').mkdir()
    geometry=generate(inputs/'geometry')
    shutil.copy2(archive,inputs/'runtime.tar.gz'); shutil.copy2(release_path,inputs/'release.json')
    shutil.copy2(release_path.parent/'runtime/farm/start_job.sh',inputs/'start_job.sh')
    logs=f'/farm_out/{getpass.getuser()}/{workflow}'
    commands=[['swif2','create',workflow]]; jobs=[]; seeds=set()
    for kind in a.types:
        for distance in a.distances_mm:
            total=a.events_per_point if a.events_per_point is not None else a.events_per_job*a.jobs_per_point
            for chunk,start in enumerate(range(0,total,a.events_per_job)):
                count=min(a.events_per_job,total-start)
                identity=f'{a.seed}:{kind}:{distance:.12g}:{chunk}'
                seed=int.from_bytes(hashlib.sha256(identity.encode()).digest()[:8],'big')%2147483646+1
                if seed in seeds: raise ValueError('Derived seed collision; choose another --seed')
                seeds.add(seed)
                name=f'bgo{kind}_d{distance:.12g}_j{chunk:06d}'
                files={p.name:p.read_text() for p in (inputs/'geometry').glob('bdx_optics__*')
                       if f'_bgo{kind}.' in p.name or p.name.endswith('__bank.txt')}
                files['run.gcard']=(inputs/'geometry'/f'bgo{kind}.gcard').read_text()
                job=dict(format='BDX_SWIF_JOB_V1',name=name,crystal_type=kind,distance_mm=distance,
                    chunk=chunk,events=count,seed=seed,particle=a.particle,momentum_GeV=a.momentum_GeV,
                    configuration=cfg,geometry=geometry['types'][kind],files=files,
                    runtime_sha256=build['archive_sha256'])
                jobfile=inputs/(name+'.json');write_json(jobfile,job)
                command=['swif2','add-job',workflow,'-name',name,'-account',a.account,
                         '-partition',a.partition,'-constraint',a.constraint,'-cores','1',
                         '-ram',a.ram,'-time',a.time,'-disk',a.disk,
                         '-stdout',logs+'/@job_name@-@attempt_id@.out','-stderr',logs+'/@job_name@-@attempt_id@.err']
                for local,source in [('runtime.tar.gz',inputs/'runtime.tar.gz'),('start_job.sh',inputs/'start_job.sh'),('job.json',jobfile)]:
                    command+=['-input',local,'file:'+str(source)]
                for local,ext in [('optical.hipo','.hipo'),('result.json','.json'),('logs.tar.gz','.logs.tar.gz')]:
                    command+=['-output','run/'+local,'file:'+str(output/'results'/(name+ext))]
                command+=['/bin/bash','start_job.sh',build['archive_sha256'],sha(jobfile)]
                commands.append(command)
                jobs.append(dict(name=name,input=str(jobfile.relative_to(output)),events=count,crystal_type=kind,
                    distance_mm=distance,seed=seed,hipo='results/'+name+'.hipo',receipt='results/'+name+'.json'))
    commands.append(['swif2','run',workflow])
    pins={str(p.relative_to(output)):sha(p) for p in inputs.rglob('*') if p.is_file()}
    plan=dict(format='BDX_SWIF_CAMPAIGN_V1',workflow=workflow,geometry_source_sha256=sha(SOURCE),
              runtime_sha256=build['archive_sha256'],configuration=cfg,jobs=jobs,commands=commands,pins=pins)
    write_json(output/'plan.json',plan)
    (output/'commands.sh').write_text('#!/usr/bin/env bash\n# Review only: submit through plan_bdx.py --submit-plan.\n'+ '\n'.join(shlex.join(c) for c in commands)+'\n')
    print(f'Prepared {len(jobs)} jobs, {sum(j["events"] for j in jobs)} events. No jobs submitted.\nPlan: {output}/plan.json')
    print(f'Submit: python3 farm/plan_bdx.py --submit-plan {output}')
    return plan


def submit(directory):
    directory=directory.resolve(); plan=json.loads((directory/'plan.json').read_text())
    for name,digest in plan['pins'].items():
        if sha(directory/name)!=digest: raise ValueError('Campaign input changed: '+name)
    if not shutil.which('swif2'): raise ValueError('swif2 not available; submit on ifarm')
    status=directory/'submission.json'
    # Exclusive marker prevents duplicate submissions, including after a partial failure.
    with status.open('x') as f: json.dump(dict(state='started',completed_commands=0),f)
    report=dict(state='submitting',completed_commands=0)
    try:
        for command in plan['commands']:
            subprocess.run(command,check=True)
            report['completed_commands']+=1;write_json(status,report)
        report['state']='submitted';write_json(status,report)
    except BaseException as error:
        report.update(state='failed',error=str(error));write_json(status,report);raise
    print('Monitor: swif2 status '+plan['workflow'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release',type=Path);p.add_argument('--output',type=Path);p.add_argument('--account');p.add_argument('--workflow')
    p.add_argument('--submit-plan',type=Path)
    p.add_argument('--types',nargs='+',type=int,choices=range(1,9),default=list(range(1,9)))
    p.add_argument('--distances-mm',nargs='+',type=float,default=[20,80,140,200])
    p.add_argument('--events-per-job',type=int,default=5);p.add_argument('--jobs-per-point',type=int,default=10)
    p.add_argument('--events-per-point',type=int,help='exact target; final job may be smaller')
    p.add_argument('--seed',type=int,default=20261125)
    p.add_argument('--particle',choices=['mu-','mu+','e-','e+','proton','pi-','pi+','gamma'],default='mu-')
    p.add_argument('--momentum-GeV',type=float,default=4)
    p.add_argument('--config',type=Path,default=ROOT/'experiments/bdx_optics/geometry_only.conf')
    p.add_argument('--set',action='append',default=[])
    p.add_argument('--partition',default='production');p.add_argument('--constraint',default='el9')
    p.add_argument('--ram',default='2GB');p.add_argument('--disk',default='2GB');p.add_argument('--time',default='30min')
    a=p.parse_args()
    if a.submit_plan: submit(a.submit_plan)
    else: prepare(a)

if __name__=='__main__':
    try: main()
    except (OSError,ValueError,RuntimeError,subprocess.SubprocessError) as error: raise SystemExit(str(error))
