"""Reproduce v0.7 Q4/T6 transient-wave diagnostic results.

Run from repository root with the package installed or PYTHONPATH=src.
Results are diagnostic comparisons, not independent experimental validation.
"""
from dataclasses import replace
from pathlib import Path
import copy
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

from thermofem2d.dynamic_diagnostics import (
    assemble_wave_system, compare_wave, integrate_wave,
    lowest_frequencies, rms_vector,
)


def write_csv(path, rows):
    if not rows:
        raise ValueError('no rows to save')
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def static_scale(system):
    dofs = np.arange(2*len(system.mesh.nodes))
    free = np.setdiff1d(dofs, list(system.fixed_u))
    u = np.zeros(len(dofs))
    u[free] = spsolve(system.Ku[free,:][:,free].tocsc(), system.traction[free])
    return rms_vector(u.reshape(-1, 2))


def end_error(t6, q4, static, *, perturbation=False):
    z = compare_wave(t6, q4, static_scale=static, perturbation=perturbation)
    return {key: float(z[key][-1]) for key in
            ('rms_difference_m','q4_rms_displacement_m',
             'relative_difference','static_normalized_difference')}


def make_plot(x,yseries,xlabel,ylabel,title,path,*,yscale=None):
    fig, ax = plt.subplots(figsize=(7.6,4.6))
    for label,values in yseries.items():
        ax.plot(x,values,marker='o',label=label)
    ax.set(xlabel=xlabel,ylabel=ylabel,title=title)
    if yscale:ax.set_yscale(yscale)
    ax.grid(alpha=.24)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path,dpi=170)
    plt.close(fig)


def run(out=Path('results/verification')):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    config=json.loads(Path('examples/coupled_plate_t6.json').read_text(encoding='utf-8'))
    material=config['material']
    final_time=1.6e-4
    steps=64
    subdivisions=[4,6,8,12,16]
    systems={}
    def system(n,method):
        key=(n,method)
        if key not in systems:
            case=copy.deepcopy(config)
            case['geometry'].update(per_side=n,radial=n)
            systems[key]=assemble_wave_system(case,method)
        return systems[key]

    # A: Pure mechanical waves, eliminating the thermal boundary step.
    wave=[]
    for n in subdivisions:
        q,t=system(n,'Q4'),system(n,'T6')
        ref=static_scale(q)
        for case,rise in [('step',None),('cosine_0.04ms',.04e-3),
                          ('cosine_0.16ms',.16e-3)]:
            profile='step' if rise is None else 'cosine'
            rw=1.6e-4 if rise is None else rise
            qr=integrate_wave(q,steps=steps,final_time=final_time,profile=profile,
                              rise_time=rw,feedback=False,thermal_case='reference')
            tr=integrate_wave(t,steps=steps,final_time=final_time,profile=profile,
                              rise_time=rw,feedback=False,thermal_case='reference')
            wave.append(dict(subdivisions=n,load_profile=case,
                        rise_time_s=(0. if rise is None else rise),steps=steps,
                        q4_static_displacement_rms_m=ref,
                        **end_error(tr,qr,ref)))
    write_csv(out/'v07_mechanical_wave_convergence.csv',wave)
    profile_labels = {'step':'Abrupt step', 'cosine_0.04ms':'Cosine ramp (0.04 ms)',
                      'cosine_0.16ms':'Cosine ramp (0.16 ms)'}
    make_plot(subdivisions,{profile_labels[name]:[100*r['relative_difference'] for r in wave
                                      if r['load_profile']==name]
                      for name in ('step','cosine_0.04ms','cosine_0.16ms')},
              'Subdivisions along each edge and radial direction',
              'T6–Q4 RMS displacement difference / Q4 RMS (%)',
              'Mechanical-only transient: load regularity and mesh refinement',
              out/'v07_mechanical_wave_convergence.png')

    # B: Separately control thermal initialization with coupling enabled.
    cases=['reference','boundary_step','preheated_equilibrium','thermal_step_only']
    init_rows=[]
    for n in (4,8,12,16):
        q,t=system(n,'Q4'),system(n,'T6')
        scale=static_scale(q)
        for scenario in cases:
            if scenario=='thermal_step_only' and n not in (4,8,16):continue
            thermal_case='boundary_step' if scenario=='thermal_step_only' else scenario
            qwork,twork=q,t
            if scenario=='thermal_step_only':
                qwork=replace(q,traction=np.zeros_like(q.traction))
                twork=replace(t,traction=np.zeros_like(t.traction))
            qr=integrate_wave(qwork,steps=steps,final_time=final_time,profile='step',
                              feedback=True,thermal_case=thermal_case)
            tr=integrate_wave(twork,steps=steps,final_time=final_time,profile='step',
                              feedback=True,thermal_case=thermal_case)
            init_rows.append(dict(subdivisions=n,scenario=scenario,
                        comparison=('load-induced increment' if scenario=='preheated_equilibrium'
                                    else 'total displacement'),
                        steps=steps,**end_error(tr,qr,scale,
                             perturbation=scenario=='preheated_equilibrium')))
    write_csv(out/'v07_thermal_initialization.csv',init_rows)
    graph_n=[4,8,16]
    initialization_labels = {'reference':'Uniform reference T',
             'boundary_step':'Abrupt thermal boundary step',
             'preheated_equilibrium':'Preheated equilibrium (increment)',
             'thermal_step_only':'Thermal step without traction'}
    make_plot(graph_n,{initialization_labels[case]:[100*next(r['relative_difference'] for r in init_rows
                                   if r['subdivisions']==n and r['scenario']==case)
                           for n in graph_n]
                       for case in ('reference','boundary_step','preheated_equilibrium',
                                    'thermal_step_only')},
              'Subdivisions along each edge and radial direction',
              'T6–Q4 RMS displacement difference / Q4 RMS (%)',
              'Thermal initialization dominates the early dynamic discrepancy',
              out/'v07_thermal_initialization.png')

    # C: Uncoupled mechanical eigenfrequencies; independent stiffness and mass.
    modes=[]
    for n in subdivisions:
        for method in ('Q4','T6'):
            vals=lowest_frequencies(system(n,method),count=6)
            for mode,frequency in enumerate(vals,1):
                modes.append(dict(subdivisions=n,element=method,mode=mode,
                                  frequency_hz=float(frequency)))
    write_csv(out/'v07_modal_frequencies.csv',modes)
    make_plot(subdivisions,{f'{method} mode {mode}':[
                    next(r['frequency_hz'] for r in modes
                         if r['subdivisions']==n and r['element']==method and r['mode']==mode)
                      for n in subdivisions]
                    for method,mode in [('Q4',1),('T6',1),('Q4',3),('T6',3)]},
              'Subdivisions along each edge and radial direction',
              'Undamped mechanical frequency (Hz)',
              'Convergence of selected undamped Q4 and T6 eigenfrequencies',
              out/'v07_modal_frequencies.png')

    # D: Temporal convergence, reporting successive-refinement changes.
    temporal=[]
    for method in ('Q4','T6'):
        s=system(16,method)
        for profile,rise in [('step',1.6e-4),('cosine_0.16ms',1.6e-4)]:
            previous=None
            for ns in (32,64,128,256):
                r=integrate_wave(s,steps=ns,final_time=final_time,
                    profile='step' if profile=='step' else 'cosine',rise_time=rise,
                    thermal_case='reference',feedback=False)
                result=dict(element=method,subdivisions=16,load_profile=profile,
                            time_steps=ns,displacement_rms_m=rms_vector(r.displacement_history[-1]),
                            delta_vs_half_steps_m=None,
                            delta_vs_half_steps_relative=None)
                if previous is not None:
                    err=rms_vector(r.displacement_history[-1]-previous.displacement_history[-1])
                    result['delta_vs_half_steps_m']=err
                    result['delta_vs_half_steps_relative']=err/result['displacement_rms_m']
                temporal.append(result)
                previous=r
    write_csv(out/'v07_time_step_sensitivity.csv',temporal)
    t_axis=[64,128,256]
    make_plot(t_axis,{f'{m} {c}':[
                      100*next(r['delta_vs_half_steps_relative'] for r in temporal
                          if r['element']==m and r['load_profile']==c and r['time_steps']==ns)
                      for ns in t_axis]
                  for m in ('Q4','T6') for c in ('step','cosine_0.16ms')},
              'Number of time steps over 0.16 ms',
              'RMS change versus half as many steps / RMS displacement (%)',
              'Temporal sensitivity at n=16: step versus smooth ramp',
              out/'v07_time_step_sensitivity.png',yscale='log')

    # E: Characteristic thermal-diffusion length versus geometric mesh spacing.
    alpha=material['conductivity']/(material['density']*material['heat_capacity'])
    depth=np.sqrt(alpha*final_time)
    typical_ligament=(config['geometry']['height']/2 -
                      (abs(config['geometry']['hole_cy'])+
                       config['geometry']['hole_height']/2))
    physics=dict(thermal_diffusivity_m2_s=alpha,observation_time_s=final_time,
                 sqrt_alpha_time_m=depth,
                 reference_outer_hole_ligament_m=typical_ligament,
                 radial_spacing_at_16_subdivisions_m=typical_ligament/16,
                 ratio_radial_spacing_to_diffusion_length_at_16=(typical_ligament/16)/depth)

    # F: Representative time histories at a shared outer-boundary vertex.
    at=[]
    for scenario,profile in [('reference','cosine'),('boundary_step','step')]:
        for n in (4,16):
            for method in ('Q4','T6'):
                s=system(n,method)
                r=integrate_wave(s,steps=128,final_time=final_time,
                       profile=profile,rise_time=final_time,thermal_case=scenario,feedback=True)
                probe=np.array([0.,config['geometry']['height']/2])
                dist,index=cKDTree(r.mesh.nodes).query(probe)
                if dist>1e-12:raise ValueError('top-center probe not in mesh')
                for i,now in enumerate(r.times):
                    at.append(dict(scenario=scenario,subdivisions=n,element=method,
                                  time_s=float(now),top_center_uy_m=float(r.displacement_history[i,index,1]),
                                  temperature_K=float(r.temperature_history[i,index])))
    write_csv(out/'v07_probe_history.csv',at)
    for scenario in ('reference','boundary_step'):
        fig,ax=plt.subplots(figsize=(7.6,4.6))
        for n in (4,16):
            for method in ('Q4','T6'):
                subset=[r for r in at if r['scenario']==scenario and
                        r['subdivisions']==n and r['element']==method]
                ax.plot([1e3*r['time_s'] for r in subset],
                        [1e6*r['top_center_uy_m'] for r in subset],label=f'{method}, n={n}')
        ax.set(xlabel='Time (ms)',ylabel='Top-center vertical displacement (µm)',
               title=('Smooth mechanical ramp, initially uniform temperature'
                      if scenario=='reference'
                      else 'Abrupt mechanical load and abrupt thermal boundary step'))
        ax.grid(alpha=.24);ax.legend();fig.tight_layout()
        fig.savefig(out/f'v07_probe_{scenario}.png',dpi=170);plt.close(fig)

    summary={
        'status':'verification study; independent external validation outstanding',
        'model':'Small-strain linear plane-stress transient thermoelasticity',
        'reference':'Independently assembled Q4 spatial FEM with matching Newmark/BE integrator',
        'mesh_divisions':subdivisions,
        'constant_top_traction_Pa':2e6,
        'final_time_s':final_time,
        'spatial_comparison_steps':steps,
        'thermal_initialization':physics,
        'limits':['No experimental or external commercial-FEM validation',
                  'Abrupt thermal boundary step is unresolved on available spatial meshes',
                  'The n=16 abrupt-load result requires additional temporal refinement',
                  'The preheated case compares incremental displacement about distinct equilibria']
    }
    (out/'v07_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({'spatial_rows':len(wave),'initialization_rows':len(init_rows),
                      'modal_rows':len(modes),'temporal_rows':len(temporal),
                      'probe_rows':len(at), 'physics':physics},indent=2))
    return summary


if __name__=='__main__':
    run()
