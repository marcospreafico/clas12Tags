#!/usr/bin/env python3
"""Export BDX HIPO banks to a new CSV directory without modifying production data."""
import argparse
from pathlib import Path
from bdx_data import export_hipo

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('hipo',type=Path)
    p.add_argument('--output',required=True,type=Path,help='new directory inside clone for Docker')
    p.add_argument('--native',action='store_true')
    a=p.parse_args()
    export_hipo(a.hipo,a.output,a.native)
    print(f'Exported Config.csv, Event.csv, Edep.csv, Flux.csv and banks.json to {a.output}')
