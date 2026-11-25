#!/usr/bin/env python3
"""Roundtrip known photon positions through HIPO export and offline area selection."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_bdx_optics import analyze
from bdx_data import IMAGE

with tempfile.TemporaryDirectory(prefix='postselection_',dir=ROOT/'build/bdx_optics') as temporary:
    root=Path(temporary); campaign=root/'production'; run=campaign/'bgo1/d120mm'; run.mkdir(parents=True)
    source=ROOT/'build/bdx_optics/transport_tests/large_event.hipo'
    target=run/'optical.hipo'; shutil.copy2(source,target)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=dict(status='completed',geometry=dict(collector='full_large_base_first_arrival'),
        configuration=dict(light_yield_per_MeV='100000',wrap_reflectivity='0.98',absorption_length_mm='10000'),
        image=IMAGE,events_per_point=2,
        runs=[dict(status='completed',crystal_type=1,distance_mm=120,directory='bgo1/d120mm',
                   hipo_file='optical.hipo',hipo_sha256=sha(target))])
    plan=campaign/'campaign.json';plan.write_text(json.dumps(manifest))
    before={str(p.relative_to(campaign)):sha(p) for p in campaign.rglob('*') if p.is_file()}
    full=analyze(campaign,root/'full')
    square=analyze(campaign,root/'square',square_mm=3)
    assert full[0]['events']==square[0]['events']==2
    assert full[0]['generated_scint']==square[0]['generated_scint']==100001
    assert full[0]['arrived_scint']==100000
    assert square[0]['arrived_scint']==50000
    assert full[0]['edep_MeV']==square[0]['edep_MeV']==2
    with (root/'square/banks/bgo1/d120mm/Event.csv').open() as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==2 and int(rows[1]['arrived_scint'])==0
    assert before=={str(p.relative_to(campaign)):sha(p) for p in campaign.rglob('*') if p.is_file()}
    result=dict(status='passed',events=2,generated=100001,full_base=100000,square_3mm=50000,
                retained_zero_arrival_events=1,production_unchanged=True)
    (ROOT/'build/bdx_optics/transport_tests/postselection.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
