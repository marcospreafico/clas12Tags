#!/usr/bin/env python3
"""Aggregate scintillation arrivals with generated-photon and MeV normalisation."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime, timezone
from bdx_data import export_hipo


def summarize(rows, light_yield):
    if not rows: raise ValueError('No events')
    n=len(rows)
    generated=sum(int(r['generated_scint']) for r in rows)
    arrived=sum(int(r['arrived_scint']) for r in rows)
    energy=sum(float(r['edep_MeV']) for r in rows)
    if generated<=0 or energy<=0: raise ValueError('Nonpositive scintillation denominator')
    efficiency=arrived/generated
    if not 0<=efficiency<=1: raise ValueError('Unphysical collection efficiency')
    # Ratio-of-totals uncertainty, treating each particle event as a cluster.
    # This includes event-to-event deposition/secondary correlations. It is not
    # the binomial error for independent photons from a common source position.
    def ratio_sem(numerator, denominator, ratio):
        if n<2: return None
        residuals=[float(r[numerator])-ratio*float(r[denominator]) for r in rows]
        total=sum(float(r[denominator]) for r in rows)
        return math.sqrt(n/(n-1)*sum(x*x for x in residuals))/total
    per_MeV=arrived/energy
    # Wilson interval documents photon counting precision, even with zero arrivals.
    center=(efficiency+0.5/generated)/(1+1/generated)
    radius=math.sqrt(efficiency*(1-efficiency)/generated+0.25/generated**2)/(1+1/generated)
    mean_distance=sum(float(r['edep_MeV'])*float(r['edep_axial_distance_mm'])
                      for r in rows if float(r['edep_MeV'])>0)/energy
    return dict(events=n,edep_MeV=energy,generated_scint=generated,arrived_scint=arrived,
        collection_efficiency=efficiency,
        efficiency_event_sem=ratio_sem('arrived_scint','generated_scint',efficiency),
        efficiency_photon_wilson68_low=max(0,center-radius),
        efficiency_photon_wilson68_high=min(1,center+radius),
        arrived_per_MeV=per_MeV,
        arrived_per_MeV_event_sem=ratio_sem('arrived_scint','edep_MeV',per_MeV),
        generated_per_MeV=generated/energy,
        arrived_per_MeV_divided_by_input_yield=per_MeV/light_yield,
        edep_axial_distance_mm=mean_distance,
        max_optical_steps=max(int(r['max_optical_steps']) for r in rows),
        boundary_repairs=sum(int(r['boundary_repairs']) for r in rows) if all('boundary_repairs' in r for r in rows) else None,
        max_boundary_shift_mm=max(float(r['max_boundary_shift_mm']) for r in rows) if all('max_boundary_shift_mm' in r for r in rows) else None,
        secondary_birth_repairs=sum(int(r['secondary_birth_repairs']) for r in rows) if all('secondary_birth_repairs' in r for r in rows) else None,
        max_birth_shift_mm=max(float(r['max_birth_shift_mm']) for r in rows) if all('max_birth_shift_mm' in r for r in rows) else None)


def analyze(campaign, output=None, square_mm=None, native=False):
    manifest=json.loads((campaign/'campaign.json').read_text())
    if manifest['status']!='completed': raise ValueError('Campaign is not completed')
    full_base=manifest['geometry'].get('collector')=='full_large_base_first_arrival'
    if square_mm is not None and (not math.isfinite(square_mm) or square_mm<=0):
        raise ValueError('Square side must be positive and finite')
    if square_mm is not None and not full_base:
        raise ValueError('Old small-SiPM productions do not contain full-base flux coordinates; run a new production')
    label='full_base' if full_base and square_mm is None else f'square_{square_mm:g}mm' if square_mm else 'legacy_sipm'
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output=(output or campaign.parent/(campaign.name+'_analysis')/(label+'_'+stamp)).resolve()
    if output==campaign or output.is_relative_to(campaign):
        raise ValueError('Analysis output must be outside the production directory')
    output.mkdir(parents=True,exist_ok=False)
    data=[]
    for run in manifest['runs']:
        if run['status']!='completed': raise ValueError('Incomplete point')
        directory=campaign/run['directory']
        if full_base:
            path=directory/run['hipo_file']
            if digest(path)!=run['hipo_sha256']: raise ValueError(f'HIPO changed after production: {path}')
            exported=output/'banks'/run['directory']
            export_hipo(path,exported,native,manifest['image'])
            with (exported/'Event.csv').open() as stream: rows=list(csv.DictReader(stream))
            if square_mm is not None:
                selected=[0]*len(rows)
                with (exported/'Flux.csv').open() as stream:
                    for photon in csv.DictReader(stream):
                        if int(photon['origin'])!=1: continue
                        if abs(float(photon['lx_mm']))<=square_mm/2 and abs(float(photon['ly_mm']))<=square_mm/2:
                            selected[int(photon['record_index'])]+=1
                # Preserve the total budget check before applying the analysis mask.
                for row,count in zip(rows,selected): row['selected_scint']=count
        else:
            path=directory/'events.csv'
            if digest(path)!=run['events_sha256']: raise ValueError(f'Event CSV changed after run: {path}')
            with path.open() as stream: rows=list(csv.DictReader(stream))
        if len(rows)!=manifest['events_per_point']: raise ValueError('Missing events')
        for row in rows:
            generated=sum(int(row[k]) for k in ('generated_scint','generated_cherenkov','injected_optical'))
            terminated=sum(int(row[k]) for k in ('arrived_scint','arrived_cherenkov','arrived_injected',
                                                 'bulk_absorbed','surface_absorbed','escaped'))
            if generated!=terminated: raise ValueError('Optical budget mismatch')
            if 'selected_scint' in row: row['arrived_scint']=row['selected_scint']
        data.append(dict(crystal_type=run['crystal_type'],distance_mm=run['distance_mm'],
            **summarize(rows,float(manifest['configuration']['light_yield_per_MeV']))))
    with (output/'collection.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    selection=dict(area='square' if square_mm is not None else 'full_base' if full_base else 'legacy_sipm',
                   square_side_mm=square_mm,center='large-base area centroid',origin='scintillation')
    plotting_version=svg_plot(data,output/'collection.svg',manifest,selection)
    (output/'analysis.json').write_text(json.dumps(dict(
        production=str(campaign),selection=selection,
        normalisation='ratio of total selected scintillation arrivals to total generated scintillation photons',
        uncertainty='event cluster SEM; photon Wilson68 interval is counting precision only',
        units=dict(distance='mm',energy='MeV'),
        analysis_code_sha256=digest(Path(__file__)),
        plotting_library='matplotlib',plotting_version=plotting_version,
        campaign_sha256=digest(campaign/'campaign.json'),points=data),indent=2)+'\n')
    print(f'Analysis: {output}')
    return data


def digest(path):
    checksum=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024),b''): checksum.update(block)
    return checksum.hexdigest()


def svg_plot(data, path, manifest, selection):
    import os
    import tempfile
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'bdx-matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#444444','#882255']
    fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
    for kind in sorted({p['crystal_type'] for p in data}):
        points=sorted((p for p in data if p['crystal_type']==kind),key=lambda p:p['distance_mm'])
        ax.errorbar([p['distance_mm'] for p in points],
                    [p['collection_efficiency'] for p in points],
                    yerr=[p['efficiency_event_sem'] or 0 for p in points],
                    marker='o',ms=4,capsize=3,color=colors[kind-1],label=f'Tipo {kind}')
    ax.set(xlabel='Distanza assiale dalla base maggiore [mm]',
           ylabel=r'$N_{\mathrm{arrivati}} / N_{\mathrm{generati}}$ (scintillazione)',
           xlim=(0,240))
    cfg=manifest['configuration']
    if selection['area']=='full_base': area='Base maggiore intera'
    elif selection['area']=='square':
        side=selection['square_side_mm']; area=f'Selezione {side:g} × {side:g} mm²'
    else:
        side=manifest['geometry']['sipm_side_mm']; area=f'SiPM simulato {side:g} × {side:g} mm²'
    label=(f'{area}; R={float(cfg["wrap_reflectivity"]):g}; '
           f'L ass.={float(cfg["absorption_length_mm"])/1000:g} m; '
           f'{manifest["events_per_point"]} eventi/punto; '
           f'{float(cfg["light_yield_per_MeV"]):g} fotoni/MeV')
    ax.set_title(label,fontsize=10,pad=10)
    fig.suptitle('BDX BGO — raccolta ottica',fontsize=14)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25)
    ax.legend(title='Geometria',bbox_to_anchor=(1.01,1),loc='upper left')
    fig.savefig(path)
    fig.savefig(path.with_suffix('.png'),dpi=170)
    plt.close(fig)
    return matplotlib.__version__


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('--output',type=Path,help='new analysis directory outside production')
    parser.add_argument('--square-mm',type=float,help='post-select a square at the large-base area centroid')
    parser.add_argument('--native',action='store_true')
    args=parser.parse_args()
    data=analyze(args.campaign.resolve(),args.output,args.square_mm,args.native)
    print(f'Analyzed {len(data)} points; collection.csv, analysis.json, collection.svg, collection.png')
