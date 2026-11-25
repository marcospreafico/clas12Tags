#!/usr/bin/env python3
"""End-to-end local staging test; never invokes SWIF2 or contacts the farm."""
import csv,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
release=Path(sys.argv[1]).resolve()
meta=json.loads(release.read_text())
base=ROOT/'build/bdx_optics/farm_test'; campaign=base/'campaign'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
if '--analyze-only' not in sys.argv[2:]:
    subprocess.run([sys.executable,str(ROOT/'farm/plan_bdx.py'),'--release',str(release),
        '--account','local_test','--workflow','bdx_local_staging','--output',str(campaign),
        '--types','1','8','--distances-mm','20','--events-per-job','2','--events-per-point','3',
        '--set','light_yield_per_MeV=80'],check=True)
    plan=json.loads((campaign/'plan.json').read_text()); sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for job in plan['jobs']:
        stage=base/'staging'/job['name'];stage.mkdir(parents=True)
        shutil.copy2(campaign/'inputs/runtime.tar.gz',stage/'runtime.tar.gz')
        shutil.copy2(campaign/'inputs/start_job.sh',stage/'start_job.sh')
        # SWIF may stage inputs as symlinks: workers must write under cwd, not the input cache.
        (stage/'job.json').symlink_to(campaign/job['input'])
        with (stage/'worker.log').open('w') as log:
            subprocess.run(['/bin/bash','start_job.sh',meta['archive_sha256'],sha(campaign/job['input'])],
                cwd=stage,stdout=log,stderr=subprocess.STDOUT,check=True)
        for local,remote in [('optical.hipo',job['hipo']),('result.json',job['receipt']),
                             ('logs.tar.gz','results/'+job['name']+'.logs.tar.gz')]:
            shutil.copy2(stage/'run'/local,campaign/remote)
        print('Staged/validated:',job['name'],flush=True)
docker_analysis='--docker-analysis' in sys.argv[2:]
reader=ROOT/'build/bdx_optics/bdx_hipo' if docker_analysis else release.parent/'runtime/build/bdx_optics/bdx_hipo'
backend=[] if docker_analysis else ['--native']
for label,extra in [('full',[]),('square',['--square-mm','3'])]:
    subprocess.run([sys.executable,str(ROOT/'scripts/analyze_bdx_hipos.py'),*backend,'--reader',str(reader),
        '--campaign',str(campaign),'--output',str(base/('analysis_'+label+('_docker' if docker_analysis else ''))),*extra],check=True)
full=json.loads((base/('analysis_full_docker/analysis.json' if docker_analysis else 'analysis_full/analysis.json')).read_text())
square=json.loads((base/('analysis_square_docker/analysis.json' if docker_analysis else 'analysis_square/analysis.json')).read_text())
assert len(full['sources'])==4 and len(full['points'])==2
assert sum(p['events'] for p in full['points'])==6
assert all(p['events']==3 and p['files']==2 for p in full['points'])
assert all(p['generated_scint']==q['generated_scint'] and p['edep_MeV']==q['edep_MeV'] and
           q['arrived_scint']<=p['arrived_scint'] for p,q in zip(full['points'],square['points']))
for entry in full['sources']:assert sha(entry['path'])==entry['sha256']
# Large photon banks: reducer must retain a zero-arrival event and exact square acceptance.
fixture=ROOT/'build/bdx_optics/transport_tests/large_event.hipo'
if fixture.exists():
    for label,extra,expected in [('full',[],100000),('square',['3'],50000)]:
        report=base/('large_'+label+'.json')
        from argparse import Namespace
        from analyze_bdx_hipos import reduce_file
        reduce_file(fixture,report,Namespace(reader=reader,native=not docker_analysis,square_mm=3 if extra else None))
        rows=json.loads(report.read_text())['events']
        assert len(rows)==2 and sum(r['selected_scint'] for r in rows)==expected
        assert rows[1]['selected_scint']==0
summary=dict(status='passed',jobs=4,events=6,points=2,symlink_input_staging=True,
    full_base_arrivals=sum(p['arrived_scint'] for p in full['points']),
    square_3mm_arrivals=sum(p['arrived_scint'] for p in square['points']),
    generated_scint=sum(p['generated_scint'] for p in full['points']),
    raw_data_unchanged=True,large_bank_reduction=fixture.exists(),swif2_submitted=False)
(base/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
