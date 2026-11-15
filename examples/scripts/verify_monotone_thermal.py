"""Reproduce v0.7 discrete-maximum-principle and heat-balance audit.

Produces traceable CSV, JSON and PNG results. No public application routes
are modified, and the mixed-order thermal reference is labelled explicitly.
"""
from pathlib import Path
import csv
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from thermofem2d.thermal_shock_reference import (
    solve_thermal_shock, shock_error, halfspace_temperature,
    extract_profile, assemble_graded_wave_system,
)
from thermofem2d.dynamic_diagnostics import integrate_wave, compare_wave


def config():
    cfg=json.loads(Path('examples/coupled_plate_t6.json').read_text(encoding='utf-8'))
    cfg['geometry'].update(width=2.,height=2.,hole_width=.24,hole_height=.20,
        hole_cx=.12,hole_cy=-.08,per_side=6,radial=28,element_type='T6')
    cfg['material'].update(conductivity=54.,density=7830.,
        heat_capacity=452.,thickness=.02,T_ref=298.)
    cfg['analysis'].update(displacement_bc={'ux':{'left':0}, 'uy':{'bottom':0}})
    return cfg


def csv_rows(path,rows):
    with path.open('w',newline='',encoding='utf-8') as fp:
        writer=csv.DictWriter(fp,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def run(out_dir='results/verification'):
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    cfg=config()
    baseline=[]
    results={}
    for fam in ('T3','Q4','T6'):
        for scheme in ('consistent','monotone'):
            res=solve_thermal_shock(cfg,fam,growth=1.45,steps=128,thermal_scheme=scheme)
            e=shock_error(res,cfg)
            baseline.append(dict(family=fam,scheme=scheme,nodes=len(res.mesh.nodes),
                elements=len(res.mesh.cells if fam=='Q4' else res.mesh.elements),
                radial=28,per_side=6,steps=128,
                first_radial_spacing_m=e['first_radial_spacing_m'],
                nodal_min_temperature_K=res.min_temperature_all_steps_K,
                nodal_max_temperature_K=res.max_temperature_all_steps_K,
                undershoot_K=e['undershoot_K'],overshoot_K=e['overshoot_K'],
                near_edge_rms_error_K=e['rms_temperature_error_K'],
                surface_flux_relative_error=e['surface_flux_relative_error'],
                boundary_heat_input_J=res.boundary_energy_input_J,
                stored_heat_change_J=res.stored_energy_change_J,
                discrete_heat_balance_abs_error_J=res.heat_balance_absolute_error_J))
            results[(fam,scheme)]=res
    csv_rows(out/'v07_monotone_thermal_comparison.csv',baseline)
    plt.figure(figsize=(8,4.7))
    pos=np.arange(3); width=.34
    for i,scheme in enumerate(('consistent','monotone')):
        sub=[x for x in baseline if x['scheme']==scheme]
        plt.bar(pos+(i-.5)*width,[x['undershoot_K'] for x in sub],width,label=scheme.title())
    plt.xticks(pos,['T3','Q4','T6-node / P1-subtriangle thermal'])
    plt.ylabel('Global temperature undershoot (K)')
    plt.title('Pure-conduction discrete maximum principle (28 graded layers)')
    plt.legend();plt.grid(axis='y',alpha=.2);plt.tight_layout()
    plt.savefig(out/'v07_monotone_undershoot_comparison.png',dpi=175);plt.close()

    depths=np.linspace(0,0.00025,251)
    exact=halfspace_temperature(depths,1.6e-4,cfg['material'])
    plt.figure(figsize=(8,4.7))
    plt.plot(depths*1000,exact,linewidth=2,label='Analytical half-space')
    for name,style in [('consistent','--'),('monotone','-.')]:
        r=results[('T6',name)]
        plt.plot(depths*1000,extract_profile(r,cfg['geometry'],depths),
                 linestyle=style,label=f'T6-node mesh: {name} thermal method')
    plt.xlabel('Depth inward from heated face (mm)')
    plt.ylabel('Temperature (K)')
    plt.title('Local thermal profile; same 28-layer graded mesh')
    plt.grid(alpha=.25);plt.legend();plt.tight_layout()
    plt.savefig(out/'v07_monotone_T6_profile.png',dpi=175);plt.close()

    steps=[]
    for nt in (32,64,128,256):
        for scheme in ('consistent','monotone'):
            r=solve_thermal_shock(cfg,'T6',growth=1.45,steps=nt,thermal_scheme=scheme)
            e=shock_error(r,cfg)
            steps.append(dict(family='T6',scheme=scheme,steps=nt,
                near_edge_rms_error_K=e['rms_temperature_error_K'],
                surface_flux_relative_error=e['surface_flux_relative_error'],
                undershoot_K=e['undershoot_K'],
                heat_balance_abs_error_J=r.heat_balance_absolute_error_J))
    csv_rows(out/'v07_monotone_timestep_sensitivity.csv',steps)
    plt.figure(figsize=(8,4.7))
    for scheme in ('consistent','monotone'):
        sub=[x for x in steps if x['scheme']==scheme]
        plt.plot([r['steps'] for r in sub],
            [r['near_edge_rms_error_K'] for r in sub],'-o',label=scheme.title())
    plt.xlabel('Time steps over 0.16 ms')
    plt.ylabel('Near-edge temperature RMS error (K)')
    plt.title('T6-node thermal transient: temporal sensitivity')
    plt.legend();plt.grid(alpha=.25);plt.tight_layout()
    plt.savefig(out/'v07_monotone_timestep_sensitivity.png',dpi=175);plt.close()

    coupled=[]
    for n,rad,nt in ((4,22,32),(6,28,64),(8,28,64)):
        cfg['geometry'].update(per_side=n,radial=rad)
        for scheme in ('consistent','monotone'):
            q=integrate_wave(assemble_graded_wave_system(cfg,'Q4',1.45,
                thermal_scheme=scheme),steps=nt,final_time=1.6e-4,
                profile='step',thermal_case='boundary_step',feedback=True)
            t=integrate_wave(assemble_graded_wave_system(cfg,'T6',1.45,
                thermal_scheme=scheme),steps=nt,final_time=1.6e-4,
                profile='step',thermal_case='boundary_step',feedback=True)
            d=float(compare_wave(t,q)['relative_difference'][-1])
            coupled.append(dict(per_side=n,radial=rad,steps=nt,scheme=scheme,
                displacement_relative_difference=d,
                Q4_min_temperature_K=float(q.temperature_history.min()),
                T6_min_temperature_K=float(t.temperature_history.min())))
    csv_rows(out/'v07_monotone_coupled_displacement.csv',coupled)
    plt.figure(figsize=(8,4.7))
    for scheme in ('consistent','monotone'):
        sub=[r for r in coupled if r['scheme']==scheme]
        plt.plot([r['per_side'] for r in sub],
                 [100*r['displacement_relative_difference'] for r in sub],'-o',
                 label=scheme.title())
    plt.xlabel('Perimeter subdivisions per side')
    plt.ylabel('T6–Q4 relative displacement difference (%)')
    plt.title('Coupled thermal shock: method sensitivity')
    plt.legend();plt.grid(alpha=.25);plt.tight_layout()
    plt.savefig(out/'v07_monotone_coupled_comparison.png',dpi=175);plt.close()
    metadata={
        'scope':'experimental numerical diagnostic; not design validated',
        'validation_status':'experimental; not independently design validated',
        'thermal_case':'pure conduction sudden 310 K left edge; 298 K right edge',
        'absolute_reference_temperature_K':298,
        'final_time_s':1.6e-4,
        'thermal_mesh':'geometric radial grading; growth=1.45, 28 layers',
        'monotone_T6_thermal':'P1 conforming four-subtriangle thermal operator on T6 nodes',
        'coupled_T6':'T6 displacement/P1-subtriangle temperature; reciprocal P2/P1 coupling',
        'limitation':'no fully independent continuum solution for coupled thermal shock',
        'data_files':['v07_monotone_thermal_comparison.csv',
                      'v07_monotone_timestep_sensitivity.csv',
                      'v07_monotone_coupled_displacement.csv'],
    }
    (out/'v07_monotone_reproducibility.json').write_text(
        json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print('Pure conduction:')
    for row in baseline:
        print(f" {row['family']:2s} {row['scheme']:10s} under={row['undershoot_K']:.4f} K, "
              f"profile_RMSE={row['near_edge_rms_error_K']:.5f} K, "
              f"flux_err={100*row['surface_flux_relative_error']:.3f}%, "
              f"heat_balance_error={row['discrete_heat_balance_abs_error_J']:.2g} J")
    print('Coupled comparisons:')
    for row in coupled:
        print(row)
    return baseline,steps,coupled

if __name__=='__main__':
    run()
