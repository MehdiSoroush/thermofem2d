"""Reproduce continuum manufactured coupled T3/T6/Q4 verification.

Outputs CSV and PNG plots; produces separate diagnostics on the previously
reported abrupt-load T6/Q4 disagreement, including a static comparison.
"""
from pathlib import Path
import copy
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

from thermofem2d.manufactured import verify_manufactured_coupled
from thermofem2d.q4_reference import (
    quad_mesh,assemble_quad,solve_q4_coupled,_edge_load,_fixed,_reduced_dirichlet,
)
from thermofem2d.geometry import make_mesh
from thermofem2d.fem import assemble_elastic, prescribed_displacements, solve_dirichlet
from thermofem2d.solve import solve_model


def write_csv(path, records):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)


def manufactured(out,base):
    rows=[]
    for n in (3,4,6,8):
        for kind in ('T3','Q4','T6'):
            r=verify_manufactured_coupled(base,kind,n,steps=8)
            rows.append({
                'element':kind,'subdivisions':n,'nodes':r.nodes,'elements':r.elements,
                'steps':r.steps,'temperature_rms_error_K':r.temperature_rms_error_K,
                'displacement_rms_error_m':r.displacement_rms_error_m,
                'temperature_max_error_K':r.temperature_max_error_K,
                'displacement_max_error_m':r.displacement_max_error_m,
                'temperature_span_exact_K':r.thermal_scale_K,
                'max_displacement_exact_m':r.displacement_scale_m,
            })
    write_csv(out/'coupled_manufactured_convergence.csv',rows)
    accelerated=[]
    for kind in ('T3','Q4','T6'):
        r=verify_manufactured_coupled(base,kind,4,steps=8,curvature=5000.)
        accelerated.append({'element':kind,'subdivisions':4,'time_curvature_per_s':5000.0,
                            'temperature_rms_error_K':r.temperature_rms_error_K,
                            'displacement_rms_error_m':r.displacement_rms_error_m})
    write_csv(out/'coupled_manufactured_acceleration.csv',accelerated)
    for key,ylabel,name in [
        ('temperature_rms_error_K','RMS nodal temperature error (K)','coupled_manufactured_temperature.png'),
        ('displacement_rms_error_m','RMS nodal displacement vector error (m)','coupled_manufactured_displacement.png'),
    ]:
        fig,ax=plt.subplots(figsize=(8,4.7))
        for kind in ('T3','Q4','T6'):
            subset=[r for r in rows if r['element']==kind]
            ax.plot([r['subdivisions'] for r in subset],[r[key] for r in subset],marker='o',label=kind)
        ax.set(yscale='log',xlabel='Subdivisions per plate side / radial direction',ylabel=ylabel)
        ax.legend();ax.grid(alpha=.25)
        fig.tight_layout();fig.savefig(out/name,dpi=160);plt.close(fig)
    return rows


def abrupt_load_diagnostics(out,base):
    """Contrast static similarity with unresolved short-time dynamics."""
    rows=[];mat=base['material']
    for n in (4,6,8,12):
        g=dict(base['geometry'],per_side=n,radial=n)
        qm=quad_mesh(g)
        _,_,Kq,_,_=assemble_quad(qm,mat)
        fq=_edge_load(qm,{'top':[0.,2e6]},mat['thickness'],2)
        fixedq={}
        for name,component in (('ux',0),('uy',1)):
            fixedq.update(_fixed(qm,base['analysis']['displacement_bc'][name],component))
        uq=_reduced_dirichlet(Kq,fq,fixedq).reshape(-1,2)
        tm=make_mesh(**dict(g,element_type='T6'))
        Kt,ft=assemble_elastic(tm,mat,np.full(len(tm.nodes),mat['T_ref']),{'top':[0.,2e6]})
        ut=solve_dirichlet(Kt,ft,prescribed_displacements(tm,base['analysis']['displacement_bc'])).reshape(-1,2)
        dist,idx=cKDTree(tm.nodes).query(qm.nodes)
        if dist.max()>1e-12:raise ValueError('mismatched static comparison nodes')
        rms=lambda v:float(np.sqrt(np.mean(np.sum(v*v,axis=1))))
        rows.append({'case':'static','subdivisions':n,'time_steps':0,
                     'relative_displacement_difference':rms(ut[idx]-uq)/rms(uq),
                     'rms_displacement_difference_m':rms(ut[idx]-uq)})
    for n in (4,6):
        for steps in (8,16,32,64):
            config=copy.deepcopy(base)
            config['geometry'].update(per_side=n,radial=n)
            config['analysis'].update(steps=steps,time_step=1.6e-4/steps)
            qref=solve_q4_coupled(config)
            config['geometry']['element_type']='T6'
            tri=solve_model(config)
            dist,idx=cKDTree(tri.mesh.nodes).query(qref.mesh.nodes)
            if dist.max()>1e-12:raise ValueError('mismatched transient comparison nodes')
            u4=qref.displacement.reshape(-1,2)
            u6=tri.displacement.reshape(-1,2)[idx]
            rows.append({'case':'abrupt_transient','subdivisions':n,'time_steps':steps,
                         'relative_displacement_difference':rms(u6-u4)/rms(u4),
                         'rms_displacement_difference_m':rms(u6-u4)})
    write_csv(out/'q4_t6_static_dynamic_diagnostics.csv',rows)
    fig,ax=plt.subplots(figsize=(8,4.7))
    static=[r for r in rows if r['case']=='static']
    ax.plot([r['subdivisions'] for r in static],
            [100*r['relative_displacement_difference'] for r in static],'-o',label='Static traction, increasing mesh density')
    ax.set(xlabel='Subdivisions per side',ylabel='RMS displacement difference relative to Q4 (%)')
    ax.grid(alpha=.25);ax.legend(loc='upper right')
    fig.tight_layout();fig.savefig(out/'q4_t6_static_comparison.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.7))
    for n in (4,6):
        sub=[r for r in rows if r['case']=='abrupt_transient' and r['subdivisions']==n]
        ax.plot([r['time_steps'] for r in sub],
                [100*r['relative_displacement_difference'] for r in sub],'-o',label=f'n={n} subdivisions')
    ax.set(xlabel='Time steps over fixed 0.16 ms window',ylabel='T6–Q4 relative RMS displacement difference (%)')
    ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(out/'q4_t6_dynamic_time_sensitivity.png',dpi=160);plt.close(fig)
    return rows


def evaluate(out=Path('results/verification')):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    base=json.loads(Path('examples/coupled_plate_t6.json').read_text())
    manufactured_rows=manufactured(out,base)
    diagnostic_rows=abrupt_load_diagnostics(out,base)
    print('T6 exact quadratic manufactured RMS temperature:',
          [r['temperature_rms_error_K'] for r in manufactured_rows if r['element']=='T6'])
    print('Abrupt-load T6/Q4 time-refinement relative differences:',
          [r['relative_displacement_difference'] for r in diagnostic_rows
           if r['case']=='abrupt_transient' and r['subdivisions']==4])
    return manufactured_rows,diagnostic_rows


if __name__=='__main__':
    evaluate()
