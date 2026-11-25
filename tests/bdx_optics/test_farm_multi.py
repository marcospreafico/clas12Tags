"""Offline SWIF2 planning, submission, worker, and pooled event statistics."""
import copy,importlib.util,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'farm')]
import plan_bdx
from release import sha,write_json
from analyze_bdx_hipos import Aggregator,discover
from argparse import Namespace


def row(seed=123,generated=100,arrived=10,energy=1):
    cfg=dict(schema_version=2,collector=1,crystal_type=1,seed=seed,wrap_reflectivity=.98,
        light_yield_per_MeV=100,absorption_length_mm=10000,refractive_index=2.15,half_length_mm=120,
        Bg_mm=24,Bm_mm=20,H_mm=25,bg_mm=12,bm_mm=10,h_mm=12,center_y_mm=1,plane_z_mm=120)
    e=dict(nominal_distance_mm=20,edep_MeV=energy,generated_scint=generated,arrived_scint=arrived,
        generated_cherenkov=0,arrived_cherenkov=0,injected_optical=0,arrived_injected=0,
        bulk_absorbed=generated-arrived,surface_absorbed=0,escaped=0,
        edep_axial_distance_mm=20,max_optical_steps=100)
    return dict(config=cfg,event=e,selected_scint=arrived,primaries=[dict(pid=13,px=4,py=0,pz=0,vx=-10,vy=0,vz=10,vt=0)])


class PooledTests(unittest.TestCase):
    def test_totals_not_file_averages_and_zero_arrivals(self):
        a=Aggregator(); a.add(Path('one.hipo'),'a',[row(123,100,10,1)])
        a.add(Path('two.hipo'),'b',[row(124,300,60,3),row(124,100,0,1)])
        data,cfg=a.results();d=data[0]
        self.assertEqual(d['events'],3);self.assertEqual(d['files'],2)
        self.assertAlmostEqual(d['collection_efficiency'],70/500)
        self.assertAlmostEqual(d['arrived_per_MeV'],14)
        self.assertGreater(d['efficiency_event_sem'],0)

    def test_repeated_content_seed_and_optical_mixture(self):
        for what in ['hash','seed','physics','geometry','beam']:
            a=Aggregator();a.add(Path('one'),'a',[row()]);r=row(124);digest='b'
            if what=='hash':digest='a'
            elif what=='seed':r['config']['seed']=123
            elif what=='physics':r['config']['absorption_length_mm']=9000
            elif what=='geometry':r['config']['Bg_mm']=30
            elif what=='beam':r['primaries'][0]['pid']=11
            with self.subTest(what=what),self.assertRaises(ValueError):a.add(Path('two'),digest,[r])

    def test_metadata_and_budget(self):
        for key in ['schema','collector','type','count']:
            r=row()
            if key=='schema':r['config']['schema_version']=1
            if key=='collector':r['config']['collector']=0
            if key=='type':r['config']['crystal_type']=0
            if key=='count':r['selected_scint']=11
            with self.subTest(key=key),self.assertRaises(ValueError):Aggregator().add(Path('a'),'a',[r])


class PlannerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name); release=self.base/'build';(release/'runtime/farm').mkdir(parents=True)
        (release/'runtime.tar.gz').write_bytes(b'test runtime')
        (release/'runtime/farm/start_job.sh').write_text('# test\n')
        self.release=release/'release.json'
        write_json(self.release,dict(archive=str(release/'runtime.tar.gz'),archive_sha256=sha(release/'runtime.tar.gz'),source=dict(source_sha256={})))
        self.args=Namespace(account='test_account',release=self.release,output=self.base/'campaign',
            workflow='bdx_test',types=[1,8],distances_mm=[20,200],events_per_job=3,jobs_per_point=5,
            events_per_point=7,seed=123,particle='mu-',momentum_GeV=4,
            config=ROOT/'experiments/bdx_optics/geometry_only.conf',set=['light_yield_per_MeV=800'],
            partition='production',constraint='el9',ram='2GB',disk='2GB',time='30min')

    def test_job_grid_and_exact_remainder_without_submission(self):
        with patch('subprocess.run') as run:
            plan=plan_bdx.prepare(self.args)
        run.assert_not_called()
        self.assertEqual(len(plan['jobs']),12)
        self.assertEqual(sum(j['events'] for j in plan['jobs']),28)
        self.assertEqual(len(set(j['seed'] for j in plan['jobs'])),12)
        self.assertEqual([j['events'] for j in plan['jobs'][:3]],[3,3,1])
        command=plan['commands'][1]
        self.assertIn('/farm_out/',command[command.index('-stdout')+1])
        self.assertIn('-input',command);self.assertIn('-output',command)
        self.assertEqual(plan['configuration']['absorption_length_mm'],'10000')
        job=json.loads((self.args.output/plan['jobs'][0]['input']).read_text())
        self.assertIn('bdx_flux',next(v for k,v in job['files'].items() if '__geometry_' in k))

    def test_submission_marker_and_input_changes(self):
        plan=plan_bdx.prepare(self.args)
        with patch('shutil.which',return_value='/fake/swif2'),patch('subprocess.run') as run:
            plan_bdx.submit(self.args.output)
            self.assertEqual(run.call_count,len(plan['commands']))
            with self.assertRaises(FileExistsError):plan_bdx.submit(self.args.output)
        self.assertEqual(json.loads((self.args.output/'submission.json').read_text())['state'],'submitted')

    def test_drift_rejected_before_submission(self):
        plan=plan_bdx.prepare(self.args)
        (self.args.output/'inputs/start_job.sh').write_text('changed')
        with patch('subprocess.run') as run,self.assertRaises(ValueError):plan_bdx.submit(self.args.output)
        run.assert_not_called()

    def test_incomplete_campaign_requires_explicit_partial(self):
        plan=plan_bdx.prepare(self.args)
        args=Namespace(campaign=[self.args.output],inputs=[],input_list=[],allow_partial=False)
        with self.assertRaises(ValueError):discover(args)
        job=plan['jobs'][0]
        (self.args.output/job['hipo']).write_bytes(b'x');(self.args.output/job['receipt']).write_text('{}')
        args.allow_partial=True
        files,expected,missing=discover(args)
        self.assertEqual(len(files),1);self.assertEqual(len(missing),11)

if __name__=='__main__':unittest.main()
