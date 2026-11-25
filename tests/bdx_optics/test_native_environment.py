"""Regression checks for interactive SSH shells and /work-style symlink paths."""
import json,os,shutil,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


class NativeEnvironmentTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name).resolve();self.physical=self.base/'physical';self.physical.mkdir()
        self.logical=self.base/'work';self.logical.symlink_to(self.physical,target_is_directory=True)
        binary=self.base/'bin';binary.mkdir();tool=binary/'geant4-config'
        tool.write_text('#!/bin/sh\necho 10.7.4\n');tool.chmod(0o755)
        self.env=dict(os.environ,PATH=str(binary)+os.pathsep+os.environ['PATH'])
        self.config=self.base/'site.env';self.config.write_text('export BDX_EXPECTED_GEANT4=10.7.4\n')

    def shell(self,body,config=None):
        return subprocess.run(['bash','--noprofile','--norc','-c',body,'test',str(ROOT/'farm/environment.sh'),
                               str(config or self.config),str(self.logical)],env=self.env,text=True,capture_output=True)

    def test_interactive_options_remain_off_and_errors_do_not_exit(self):
        p=self.shell('''set +e +u; set +o pipefail
export BDX_RUNTIME_ROOT="$3"
source "$1" "$2"
case $- in *e*|*u*) exit 41 ;; esac
false
: "$BDX_COMPLETION_UNSET"
printf 'ALIVE:%s\\n' "$GEMC"
''')
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn('ALIVE:'+str(self.physical),p.stdout)

    def test_batch_strict_options_remain_enabled(self):
        p=self.shell('''set -euo pipefail
export BDX_RUNTIME_ROOT="$3"
source "$1" "$2"
case $- in *e*) : ;; *) exit 41 ;; esac
case $- in *u*) : ;; *) exit 42 ;; esac
echo LOADED
false | true
echo UNREACHABLE
''')
        self.assertEqual(p.returncode,1,p.stderr)
        self.assertIn('LOADED',p.stdout);self.assertNotIn('UNREACHABLE',p.stdout)

    def test_failed_interactive_load_returns_error_without_exiting(self):
        p=self.shell('''set +e +u; set +o pipefail
export BDX_RUNTIME_ROOT="$3"
source "$1" "$2"
printf 'RETURNED:%s\\n' "$?"
echo ALIVE
''',self.base/'missing.env')
        self.assertEqual(p.returncode,0,p.stderr);self.assertIn('RETURNED:1',p.stdout);self.assertIn('ALIVE',p.stdout)

    def test_native_helper_uses_physical_release_path(self):
        farm=self.physical/'farm';farm.mkdir()
        for name in ['test_bdx.sh','environment.sh']:shutil.copy2(ROOT/'farm'/name,farm/name)
        release=farm/'build/20261125';runtime=release/'runtime'
        (runtime/'farm').mkdir(parents=True);(runtime/'scripts').mkdir()
        (release/'release.json').write_text('{}')
        shutil.copy2(self.config,runtime/'farm/ifarm.env')
        (runtime/'scripts/run_bdx_optics.py').write_text('import json,os,sys\nprint(json.dumps(dict(script=__file__,root=os.environ["BDX_RUNTIME_ROOT"],args=sys.argv[1:])))\n')
        environment=dict(self.env);environment.pop('BDX_RELEASE_JSON',None)
        p=subprocess.run(['bash',str(self.logical/'farm/test_bdx.sh'),'--events','1'],env=environment,text=True,capture_output=True)
        self.assertEqual(p.returncode,0,p.stderr)
        result=json.loads(p.stdout.splitlines()[-1])
        self.assertEqual(result['root'],str(runtime));self.assertEqual(result['script'],str(runtime/'scripts/run_bdx_optics.py'))
        self.assertIn('--native',result['args']);self.assertEqual(result['args'][-2:],['--events','1'])

    def test_runner_dry_run_invoked_through_symlink(self):
        tracked=subprocess.check_output(['git','-C',str(ROOT),'ls-files','scripts','source',
            'experiments/bdx/geometry_crs60.pl','experiments/bdx_optics/geometry_only.conf'],text=True).splitlines()
        for name in tracked:
            if Path(name).suffix in {'.a','.o'}:continue
            src=ROOT/name
            if not src.is_file():continue
            dst=self.physical/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
        subprocess.run(['git','init','-q',str(self.physical)],check=True)
        subprocess.run(['git','-C',str(self.physical),'-c','user.name=Test','-c','user.email=test@example.invalid',
            '-c','commit.gpgsign=false','-c','core.hooksPath=/dev/null','commit','-q','--allow-empty','-m','fixture'],check=True)
        p=subprocess.run([sys.executable,str(self.logical/'scripts/run_bdx_optics.py'),'--dry-run','--types','1',
            '--distances-mm','20','--events','1','--output',str(self.logical/'build/test')],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        manifest=json.loads((self.physical/'build/test/campaign.json').read_text())
        self.assertIn('scripts/run_bdx_optics.py',manifest['source_sha256'])
        self.assertEqual(len(manifest['runs']),1)


if __name__=='__main__':unittest.main()
