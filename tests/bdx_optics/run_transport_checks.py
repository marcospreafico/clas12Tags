"""Analytic ray, TIR, reflectivity, absorption and >100-step Geant4 checks."""
import csv
import json
import math
import os
import re
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from run_bdx_optics import config_values

output=ROOT/'build/bdx_optics/transport_tests'
results=[]


def check(name, mode, settings, n, expected, steps=False, first=None):
    config=output/f'{name}.conf'
    values=config_values(ROOT/'experiments/bdx_optics/geometry_only.conf',['absorption_length_mm=0',*settings])
    config.write_text(''.join(f'{k}={v}\n' for k,v in values.items()))
    path=output/f'{name}.csv'
    env=dict(os.environ,BDX_OPTICS_CONFIG=str(config),BDX_OPTICS_CSV=str(path))
    env.pop('BDX_OPTICS_ARRIVALS',None)
    with (output/f'{name}.log').open('w') as log:
        subprocess.run([str(output/'transport_test'),mode,str(n)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    with path.open() as stream: rows=list(csv.DictReader(stream))
    assert len(rows)==n
    hits=sum(int(r['arrived_injected']) for r in rows)
    for r in rows:
        assert int(r['injected_optical'])==1
        assert sum(int(r[k]) for k in ['arrived_injected','bulk_absorbed','surface_absorbed','escaped'])==1
    maximum=max(int(r['max_optical_steps']) for r in rows)
    if expected is None: pass
    elif expected in (0,1): assert hits==round(n*expected),(name,hits,n,expected)
    else: assert abs(hits-n*expected)<6*math.sqrt(n*expected*(1-expected)),(name,hits,n,expected)
    if steps: assert maximum>100,maximum
    first_counts=None
    if first is not None:
        matched=re.search(r'FIRST_STATS (\d+) (\d+) (\d+)',(output/f'{name}.log').read_text())
        assert matched
        first_counts=list(map(int,matched.groups()))
        for count, probability in zip(first_counts,first):
            tolerance=6*math.sqrt(n*probability*(1-probability))
            assert abs(count-n*probability)<=tolerance,(name,count,probability)
    result=dict(name=name,events=n,arrivals=hits,expected_fraction=expected,max_steps=maximum,status='passed',first_outcomes=first_counts)
    results.append(result)
    print(json.dumps(result),flush=True)


check('ideal_sensor','center',['wrap_model=contact','wrap_reflectivity=0'],100,1)
check('outside_aperture','outside',['wrap_model=contact','wrap_reflectivity=0'],100,0)
check('specular_one_bounce','one_reflection',['wrap_model=contact','wrap_reflectivity=0.85'],6000,0.85)
check('bare_escape','outside',['wrap_model=bare','refractive_index=1'],100,0)
check('bulk_absorption','center',['wrap_model=contact','absorption_length_mm=240'],6000,math.exp(-0.5))
check('bulk_absorption_10m','center',['wrap_model=contact','absorption_length_mm=10000'],6000,math.exp(-120/10000))
check('total_internal_reflection','one_reflection',['wrap_model=bare','wrap_reflectivity=0'],100,1)
check('bare_normal_escape','outside',['wrap_model=bare','wrap_reflectivity=0'],100,0)
check('paint_zero_R','center',['wrap_model=paint','wrap_reflectivity=0'],100,1)
check('paint_zero_R_outside','outside',['wrap_model=paint','wrap_reflectivity=0'],100,0)
check('long_ray_no_cut','many_reflections',['wrap_model=contact','wrap_reflectivity=1','absorption_length_mm=1e12'],20,1,True)
check('paint_diffuse_zone','outside',[],6000,None,first=(0.98,0,0.02))
check('paint_mixed_zone','mixed',[],6000,None,first=(0.49,0.49,0.02))
check('paint_specular_zone','far_side',[],6000,None,first=(0,0.98,0.02))
check('boundary_tolerance_small','mixed',['boundary_push_mm=1e-7'],6000,None,first=(0.49,0.49,0.02))
check('boundary_tolerance_large','mixed',['boundary_push_mm=1e-5'],6000,None,first=(0.49,0.49,0.02))
comparison=[r for r in results if r['name'] in {'paint_mixed_zone','boundary_tolerance_small','boundary_tolerance_large'}]
pool=sum(r['arrivals'] for r in comparison)/sum(r['events'] for r in comparison)
for r in comparison:
    assert abs(r['arrivals']/r['events']-pool) < 6*math.sqrt(pool*(1-pool)/r['events'])
# Force one rounded-outside secondary birth per charged-particle event. This
# exercises the real secondary preparation and Geant4 tracking, not just the
# geometric point helper.
values=config_values(ROOT/'experiments/bdx_optics/geometry_only.conf',['light_yield_per_MeV=50'])
config=output/'secondary_birth.conf'
config.write_text(''.join(f'{k}={v}\n' for k,v in values.items()))
path=output/'secondary_birth.csv'
env=dict(os.environ,BDX_OPTICS_CONFIG=str(config),BDX_OPTICS_CSV=str(path))
env.pop('BDX_OPTICS_ARRIVALS',None)
with (output/'secondary_birth.log').open('w') as log:
    subprocess.run([str(output/'transport_test'),'secondary_birth','10'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
with path.open() as stream: rows=list(csv.DictReader(stream))
assert len(rows)==10
assert 'BIRTH_INJECTIONS 10' in (output/'secondary_birth.log').read_text()
for row in rows:
    generated=int(row['generated_scint'])
    assert generated>0 and int(row['injected_optical'])==0
    assert generated==sum(int(row[k]) for k in ['arrived_scint','bulk_absorbed','surface_absorbed','escaped'])
    assert int(row['secondary_birth_repairs'])>=1
    assert 1.28e-6<float(row['max_birth_shift_mm'])<1e-4
result=dict(name='secondary_birth_outside_BGO',events=10,
            generated_scint=sum(int(r['generated_scint']) for r in rows),
            max_birth_shift_mm=max(float(r['max_birth_shift_mm']) for r in rows),status='passed')
results.append(result);print(json.dumps(result),flush=True)
(output/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
