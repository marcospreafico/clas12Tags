"""HIPO I/O process wrapper, shared by production QA and offline exports."""
import os
from pathlib import Path
import subprocess
from bdx_geometry import ROOT

IMAGE='jeffersonlab/gemc:dev-almalinux94@sha256:4dd5fa8c6714a7f9241106c07edb32f6c3e2e645ec55958b508b8f0d5f95ac74'


def hipo_command(path, mode, destination, native=False, image=None):
    executable=ROOT/'build/bdx_optics/bdx_hipo'
    if not executable.is_file(): raise ValueError('Compile the HIPO reader using scripts/build_bdx_optics.sh')
    path, destination=Path(path).resolve(),Path(destination).resolve()
    if native: return [str(executable),str(path),mode,str(destination)]
    def inside(p):
        if not p.is_relative_to(ROOT): raise ValueError('Docker HIPO paths must be inside the clone')
        return '/work/clas12Tags/'+str(p.relative_to(ROOT))
    return ['docker','run','--rm','--platform','linux/amd64','--entrypoint','/bin/bash',
            '-v',f'{ROOT}:/work/clas12Tags',image or os.environ.get('BDX_GEMC_IMAGE',IMAGE),'-lc',
            'module load hipo/4.5.0; exec "$@"','bdx-hipo',inside(executable),inside(path),mode,inside(destination)]


def export_hipo(path, destination, native=False, image=None):
    destination=Path(destination)
    command=hipo_command(path,'--export',destination,native,image)
    destination.mkdir(parents=True,exist_ok=False)
    with (destination/'export.log').open('w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
    return destination
