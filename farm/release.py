#!/usr/bin/env python3
"""Create immutable native build releases, including the worker and runtime inputs."""
import hashlib,json,os,shutil,subprocess,sys,tarfile
from datetime import datetime,timezone
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def write_json(path,value):
    path=Path(path); temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n'); temporary.replace(path)


def prepare(root,env):
    root=Path(root).resolve(); env=Path(env).resolve()
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    release=root/'farm/build'/stamp; runtime=release/'runtime'; runtime.mkdir(parents=True)
    ignore=shutil.ignore_patterns('*.o','*.a','*.hipo','gemc','.sconsign*','__pycache__','build','bin','ifarm.env')
    for name in ['source','api','scripts','farm']:
        shutil.copytree(root/name,runtime/name,ignore=ignore,symlinks=False)
    destination=runtime/'experiments/bdx'; destination.mkdir(parents=True)
    shutil.copy2(root/'experiments/bdx/geometry_crs60.pl',destination)
    destination=runtime/'experiments/bdx_optics'; destination.mkdir(parents=True)
    for p in (root/'experiments/bdx_optics').glob('*.conf'): shutil.copy2(p,destination/p.name)
    shutil.copy2(env,runtime/'farm/ifarm.env')
    git_root=subprocess.run(['git','-C',str(root),'rev-parse','--show-toplevel'],capture_output=True,text=True)
    own_git=git_root.returncode==0 and Path(git_root.stdout.strip()).resolve()==root
    commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip() if own_git else None
    status=subprocess.check_output(['git','-C',str(root),'status','--short'],text=True) if own_git else 'No git metadata; source SHA256 snapshot is authoritative.'
    write_json(release/'source.json',dict(source_root=str(root),git_commit=commit,git_status=status,
        source_sha256={str(p.relative_to(runtime)):sha(p) for p in runtime.rglob('*') if p.is_file()}))
    print(release)


def finish(release):
    release=Path(release).resolve(); runtime=release/'runtime'
    for name in ['source/gemc','build/bdx_optics/bdx_hipo']:
        if not os.access(runtime/name,os.X_OK): raise ValueError('Missing executable '+name)
    archive=release/'runtime.tar.gz'
    if archive.exists(): raise ValueError('Release already finalized')
    def keep(info):
        p=Path(info.name)
        if p.suffix in {'.o','.a','.pyc'} or p.name.startswith('.sconsign') or '__pycache__' in p.parts: return None
        return info
    with tarfile.open(archive,'w:gz') as tar: tar.add(runtime,arcname='runtime',filter=keep)
    manifest=dict(format='BDX_NATIVE_RELEASE_V1',archive=str(archive),archive_sha256=sha(archive),
        binary_sha256=sha(runtime/'source/gemc'),reader_sha256=sha(runtime/'build/bdx_optics/bdx_hipo'),
        geant4=subprocess.check_output(['geant4-config','--version'],text=True).strip(),
        modules=subprocess.run(['bash','-c','module -t list'],text=True,capture_output=True).stderr,
        source=json.loads((release/'source.json').read_text()))
    write_json(release/'release.json',manifest); print(release/'release.json')

if __name__=='__main__':
    if sys.argv[1]=='prepare': prepare(*sys.argv[2:])
    elif sys.argv[1]=='finish': finish(sys.argv[2])
    else: raise SystemExit('Use prepare ROOT ENV or finish RELEASE')
