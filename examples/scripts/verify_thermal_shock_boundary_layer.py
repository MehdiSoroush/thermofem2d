"""Reproduce the 2-D thermal-shock spatial-resolution investigation.

This script prepares a publication-review dataset, not a release commit.
The exact erfc temperature field is a *local* benchmark near the left edge.
It is not an exact solution for the full perforated-plate geometry.
"""
from pathlib import Path
from copy import deepcopy
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from thermofem2d.thermal_shock_reference import (
    halfspace_temperature, halfspace_surface_flux, thermal_diffusivity,
    penetration_length, radial_fractions, solve_thermal_shock, shock_error,
    extract_profile, assemble_graded_wave_system,
)
from thermofem2d.dynamic_diagnostics import integrate_wave, compare_wave


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results'/'verification'
T_FINAL=1.6e-4
GROWTH=1.45


def case(side=6, radial=28):
    c=json.loads((ROOT/'examples/coupled_plate_t6.json').read_text(encoding='utf-8'))
    c['geometry'].update(per_side=side,radial=radial)
    return c


def write_csv(path, rows):
    if not rows:
        raise ValueError('empty results')
    with path.open('w',newline='',encoding='utf-8') as fd:
        w=csv.DictWriter(fd,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def thermal_study():
    c=case();depth=np.linspace(0,.00025,101)
    exact=halfspace_temperature(depth,T_FINAL,c['material'])
    profile={'depth_mm': depth*1000.,'exact_halfspace_K':exact}
    rows=[];res_cache={}
    for fam in ('T3','Q4','T6'):
        for name,nr,growth in (('uniform',12,1.),('uniform28',28,1.),('graded',28,GROWTH)):
            c1=case(6,nr)
            r=solve_thermal_shock(c1,fam,growth=growth,steps=128)
            e=shock_error(r,c1)
            e['mesh_style']=name
            rows.append(e)
            res_cache[(fam,name)]=r
            profile[f'{fam}_{name}_K']=extract_profile(r,c1['geometry'],depth)
    write_csv(OUT/'v07_thermal_shock_summary.csv',rows)
    write_csv(OUT/'v07_thermal_shock_profiles.csv',[
        {name:float(val[i]) for name,val in profile.items()} for i in range(len(depth))])

    fig,ax=plt.subplots(figsize=(8.8,5.3))
    ax.plot(depth*1000.,exact,label='Semi-infinite analytical reference',lw=2.5)
    for fam,st in (('T6','uniform'),('T3','graded'),('Q4','graded'),('T6','graded')):
        ax.plot(depth*1000.,profile[f'{fam}_{st}_K'],label=f'{fam} — {st}',lw=1.6)
    ax.set(xlabel='Penetration from hot left edge (mm)',ylabel='Temperature (K)',
           title='Thermal-shock boundary layer at t = 0.16 ms')
    ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper right')
    fig.tight_layout();fig.savefig(OUT/'v07_analytical_temperature_profiles.png',dpi=190);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.8,4.5))
    for fam in ('T3','Q4','T6'):
        err=profile[f'{fam}_graded_K']-exact
        ax.plot(depth*1000,err,label=f'{fam} — graded')
    ax.axhline(0,color='0.35',lw=.8)
    ax.set(xlabel='Penetration from heated face (mm)',ylabel='Temperature difference (K)',
           title='Graded 2-D FEM minus analytical half-space reference')
    ax.grid(alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(OUT/'v07_graded_temperature_error.png',dpi=190);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.8,4.6))
    groups=['T3','Q4','T6'];x=np.arange(3);width=.34
    for name,offset in [('uniform28',-width/2),('graded',width/2)]:
        arr=[next(e['rms_temperature_error_K'] for e in rows
                  if e['family']==fam and e['mesh_style']==name) for fam in groups]
        ax.bar(x+offset,arr,width,label=('Uniform (28 layers)' if name=='uniform28' else 'Graded (28 layers)'))
    ax.set(yscale='log',xticks=x,xticklabels=groups,ylabel='RMS temperature error (K)',
           title='Local temperature accuracy: matched 28-layer meshes')
    ax.grid(axis='y',alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(OUT/'v07_thermal_resolution_summary.png',dpi=190);plt.close(fig)

    # Independently check temporal and spatial effects without assuming
    # monotonically decreasing error at each refinement level.
    temporal=[];c1=case(6,28)
    for steps in (32,64,128,256):
        r=solve_thermal_shock(c1,'T6',growth=GROWTH,steps=steps)
        temporal.append(shock_error(r,c1))
    write_csv(OUT/'v07_thermal_temporal_refinement.csv',temporal)
    spatial=[]
    for family in ('T3','Q4','T6'):
        for radial in (22,28,32):
            ci=case(6,radial)
            r=solve_thermal_shock(ci,family,growth=GROWTH,steps=128)
            spatial.append(shock_error(r,ci))
    write_csv(OUT/'v07_thermal_spatial_refinement.csv',spatial)
    return rows,profile,res_cache,temporal,spatial


def coupled_study():
    comparison=[]
    for side in (4,6,8):
        c=case(side,28)
        for name,growth in (('uniform',1.),('graded',GROWTH)):
            q=integrate_wave(
                assemble_graded_wave_system(c,'Q4',growth),
                steps=64,final_time=T_FINAL,profile='step',
                thermal_case='boundary_step',feedback=True)
            t=integrate_wave(
                assemble_graded_wave_system(c,'T6',growth),
                steps=64,final_time=T_FINAL,profile='step',
                thermal_case='boundary_step',feedback=True)
            comp=compare_wave(t,q)
            comparison.append({
                'subdivisions':side,'radial':28,'mesh_style':name,
                'growth':growth,'steps':64,'final_time_s':T_FINAL,
                'relative_displacement_difference':float(comp['relative_difference'][-1]),
                'absolute_RMS_displacement_difference_m':float(comp['rms_difference_m'][-1]),
                'Q4_RMS_displacement_m':float(comp['q4_rms_displacement_m'][-1]),
                'Q4_temperature_min_K':float(np.min(q.temperature_history[-1])),
                'T6_temperature_min_K':float(np.min(t.temperature_history[-1])),
            })
    write_csv(OUT/'v07_graded_coupled_comparison.csv',comparison)
    fig,ax=plt.subplots(figsize=(8.8,4.9))
    for name in ('uniform','graded'):
        selected=[d for d in comparison if d['mesh_style']==name]
        ax.plot([d['subdivisions'] for d in selected],
                [d['relative_displacement_difference']*100 for d in selected],
                'o-',label=f'{name.capitalize()} radial spacing')
    ax.set(xlabel='Subdivisions on each side',ylabel='RMS Q4–T6 displacement difference (%)',
           title='Coupled thermal shock: effect of boundary-layer resolution')
    ax.set_xticks((4,6,8));ax.grid(alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(OUT/'v07_coupled_shock_mesh_grading.png',dpi=190);plt.close(fig)
    return comparison


def create_heat_map(result):
    # Plot quadratic T6 temperature by subdividing the parent triangle into
    # four straight-sided subtriangles. Every midside degree of freedom is used.
    tri6=result.mesh.elements
    a,b,c,d,e,f=tri6.T
    child=np.concatenate([
        np.column_stack((a,d,f)),
        np.column_stack((d,b,e)),
        np.column_stack((f,e,c)),
        np.column_stack((d,e,f))],axis=0)
    xy=result.mesh.nodes
    x=(xy[:,0]+1)*1000
    y=xy[:,1]
    from matplotlib.tri import Triangulation
    triang=Triangulation(x,y,triangles=child)
    fig,axs=plt.subplots(1,2,figsize=(11,4.2))
    heat=axs[0].tripcolor(triang,result.temperature,shading='gouraud',vmin=298,vmax=310)
    axs[0].set(xlabel='Distance from left edge (mm)',ylabel='Vertical coordinate y (m)',
               title='T6 boundary-layer temperature')
    axs[0].set_xlim(0,.35);axs[0].set_ylim(-.35,.35)
    axs[0].set_aspect('auto')
    fig.colorbar(heat,ax=axs[0],label='Temperature (K)')
    axs[1].triplot(triang,color='0.3',lw=.35)
    axs[1].set(xlabel='Distance from left edge (mm)',
               ylabel='Vertical coordinate y (m)',title='Strong radial grading near the hot surface')
    axs[1].set_xlim(0,.35);axs[1].set_ylim(-.35,.35)
    axs[1].set_aspect('auto')
    fig.text(.5,.005,'Horizontal geometry deliberately magnified; plots are not to scale.',
             ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.025,1,1))
    fig.savefig(OUT/'v07_thermal_boundary_layer_heatmap.png',dpi=190)
    plt.close(fig)


def evaluate():
    OUT.mkdir(parents=True,exist_ok=True)
    rows,profile,cache,temporal,spatial=thermal_study()
    comparison=coupled_study()
    create_heat_map(cache[('T6','graded')])
    mat=case()['material']
    summary={
        'status':'experimental numerical verification; not design validated',
        'author':'Mehdi Soroush',
        'method':'locally graded perforated plate Q4, T3 and T6; finite heat-only FE plus coupled waves',
        'analytical_reference':'semi-infinite erfc, valid locally near left face at early time',
        'thermal_diffusivity_m2_per_s':thermal_diffusivity(mat),
        'diffusion_length_m':float(np.sqrt(thermal_diffusivity(mat)*T_FINAL)),
        'penetration_1_percent_m':float(penetration_length(mat,T_FINAL)),
        'analytical_surface_heat_flux_W_per_m2':halfspace_surface_flux(T_FINAL,mat),
        'thermal_study':rows,
        'temporal_refinement':temporal,
        'spatial_refinement':spatial,
        'coupled_study':comparison,
        'limitations':['full 2-D field has nonphysical undershoots on strongly graded meshes',
                       'semi-infinite temperature solution is a local early-time approximation',
                       'Q4/T6 coupled comparison shares time integration',
                       'no independent external FEM package or coupled continuum solution verified'],
    }
    (OUT/'v07_shock_benchmark_summary.json').write_text(
        json.dumps(summary,indent=2),encoding='utf-8')
    print('Thermal reference alpha:',summary['thermal_diffusivity_m2_per_s'])
    print('1% penetration depth in mm:',summary['penetration_1_percent_m']*1e3)
    for r in rows:
        print(r['family'],r['mesh_style'],'RMS K:',round(r['rms_temperature_error_K'],6),
              'flux rel:',round(r['surface_flux_relative_error'],4),
              'undershoot K:',round(r['undershoot_K'],4))
    for r in comparison:
        print('coupled',r['subdivisions'],r['mesh_style'],
              'relative percent:',round(r['relative_displacement_difference']*100,3))
    return summary

if __name__=='__main__':
    evaluate()
