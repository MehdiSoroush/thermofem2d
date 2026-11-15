"""Reproduce the v0.7 coupled energy balance and cooling sign study.

This numerical verification is not a publication release. The semi-discrete
operator and mixed time integrator are inherited from the existing
research harness. No public GitHub modifications are made.
"""
from pathlib import Path
import csv
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve

from thermofem2d.thermal_shock_reference import assemble_graded_wave_system
from thermofem2d.dynamic_diagnostics import integrate_wave
from thermofem2d.coupled_energy_audit import audit_coupled_energy


def write_rows(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def base_config(root):
    cfg=json.loads((root/'examples/coupled_plate_t6.json').read_text(encoding='utf-8'))
    cfg['geometry'].update(width=2.,height=2.,hole_width=.24,hole_height=.20,
        hole_cx=.12,hole_cy=-.08,per_side=4,radial=22,element_type='T6')
    cfg['material'].update(conductivity=54.,density=7830.,heat_capacity=452.,
        thickness=.02,T_ref=298.)
    cfg['analysis']['displacement_bc']={'ux':{'left':0.},'uy':{'bottom':0.}}
    return cfg


def analyze(root, out_dir=None):
    root=Path(root)
    out=Path(out_dir) if out_dir else root/'results'/'verification'
    out.mkdir(parents=True,exist_ok=True)
    cfg=base_config(root)
    records=[];saved=None
    for family in ('Q4','T6'):
        for scheme in ('consistent','monotone'):
            sys=assemble_graded_wave_system(cfg,family,1.45,thermal_scheme=scheme)
            for steps in (64,256,1024,2048):
                r=integrate_wave(sys,steps=steps,final_time=1.6e-4,profile='step',
                    thermal_case='boundary_step',feedback=True)
                eb=audit_coupled_energy(sys,r,thermal_case='boundary_step')
                records.append({
                    'element':family,'thermal_scheme':scheme,
                    'per_side':4,'radial_layers':22,'steps':steps,
                    'time_step_s':float(r.times[1]-r.times[0]),
                    'min_temperature_K':eb.min_temperature_K,
                    'min_temperature_time_s':eb.min_temperature_time_s,
                    'min_temperature_node':eb.min_temperature_node,
                    'cooling_nodes_at_final_time':eb.cooling_node_count,
                    'energy_change_J':float(eb.energy_change_J[-1]),
                    'mechanical_work_J':float(eb.mechanical_work_J[-1]),
                    'thermal_dirichlet_perturbation_work_J':float(eb.temperature_boundary_work_J[-1]),
                    'quadratic_conduction_dissipation_J':float(eb.conductivity_dissipation_J[-1]),
                    'backward_euler_numerical_dissipation_J':float(eb.backward_euler_dissipation_J[-1]),
                    'mixed_timestep_coupling_correction_J':float(eb.coupling_time_quadrature_J[-1]),
                    'max_energy_identity_residual_J':float(np.max(abs(eb.identity_residual_J))),
                    'relative_energy_identity_residual':eb.relative_identity_residual,
                    'max_free_mechanical_equation_residual_N':eb.max_mechanical_free_residual_N,
                    'max_free_thermal_equation_residual_W':eb.max_thermal_free_residual_W,
                })
                if family=='T6' and scheme=='monotone' and steps==1024:
                    saved=(sys,r,eb)
    write_rows(out/'v07_coupled_energy_timestep.csv',records)
    # Exact *material* limit for an adiabatic affine tensile strain, separate
    # from the non-adiabatic transient plate simulation.
    sign=[]
    mat=cfg['material'];eps=1.e-5
    beta=mat['E']*mat['alpha']/(1-mat['nu'])
    exact=-mat['T_ref']*beta*2*eps/(mat['density']*mat['heat_capacity'])
    for family in ('T3','Q4','T6'):
        for scheme in ('consistent','monotone'):
            sys=assemble_graded_wave_system(cfg,family,1.45,thermal_scheme=scheme)
            u=eps*sys.mesh.nodes.reshape(-1)
            nodal=-sys.reference_temperature*spsolve(sys.Ct.tocsc(),sys.G.T@u)
            sign.append({'family':family,'thermal_scheme':scheme,
                'uniform_isotropic_strain':eps,'analytic_adiabatic_cooling_K':exact,
                'nodal_temperature_change_min_K':float(nodal.min()),
                'nodal_temperature_change_max_K':float(nodal.max()),
                'max_error_K':float(max(abs(nodal-exact)))})
    write_rows(out/'v07_adiabatic_sign_check.csv',sign)
    fig,ax=plt.subplots(figsize=(8,4.7))
    for family in ('Q4','T6'):
        for scheme in ('consistent','monotone'):
            rr=[z for z in records if z['element']==family and z['thermal_scheme']==scheme]
            ax.plot([x['steps'] for x in rr],[298-x['min_temperature_K'] for x in rr],
                marker='o',label=f'{family} {scheme}')
    ax.set(xlabel='Time steps over 0.16 ms',ylabel='Minimum nodal temperature below 298 K (K)',
        title='Apparent cooling under temporal refinement (fixed graded mesh)')
    ax.set_xscale('log',base=2);ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/'v07_coupled_cooling_time_sensitivity.png',dpi=175);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.7))
    for family in ('Q4','T6'):
        for scheme in ('consistent','monotone'):
            rr=[z for z in records if z['element']==family and z['thermal_scheme']==scheme]
            ax.plot([x['steps'] for x in rr],[abs(x['mixed_timestep_coupling_correction_J']) for x in rr],
                marker='o',label=f'{family} {scheme}')
    ax.set(xlabel='Time steps over 0.16 ms',
        ylabel='Magnitude of coupled time-quadrature correction (J)',
        title='Newmark / backward-Euler coupling integration discrepancy')
    ax.set_xscale('log',base=2);ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/'v07_coupling_quadrature_time_sensitivity.png',dpi=175);plt.close(fig)
    sys,r,eb=saved
    fig,ax=plt.subplots(figsize=(8,4.7))
    scale=1e3
    ax.plot(r.times*scale,eb.energy_change_J,label='Quadratic energy change')
    ax.plot(r.times*scale,eb.mechanical_work_J+eb.temperature_boundary_work_J-
        eb.conductivity_dissipation_J-eb.backward_euler_dissipation_J+
        eb.coupling_time_quadrature_J,linestyle='--',label='Work minus losses plus time-coupling correction')
    ax.set(xlabel='Time (ms)',ylabel='Cumulative energy (J)',
        title='Independent discrete energy ledger, T6-node monotone, 1,024 steps')
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/'v07_coupled_energy_identity.png',dpi=175);plt.close(fig)
    metadata={
        'status':'experimental coupled-energy diagnostic; not continuum validated',
        'base_version':'v0.6.0 production solver with v0.7 verification studies',
        'test_problem':'2m x 2m perforated plate; 2MPa abrupt traction; left T=310K and right T=298K; T0=298K',
        'grading_factor':1.45,'final_time_s':1.6e-4,
        'energy_definition':'quadratic perturbation functional; not absolute sensible heat',
        'time_method':'Newmark average acceleration + backward Euler thermal',
        'energy_identity_note':'includes temperature-boundary work, BE dissipation, and finite-step coupling correction',
        'initial_boundary_note':'temperature step is imposed in the initial discrete state at t=0+; instantaneous heating event is outside the audited interval',
        'limitation':'time/space convergence of true coupled thermal shock not established',
        'files':['v07_coupled_energy_timestep.csv','v07_adiabatic_sign_check.csv',
            'v07_coupled_cooling_time_sensitivity.png',
            'v07_coupling_quadrature_time_sensitivity.png','v07_coupled_energy_identity.png']}
    (out/'v07_coupled_energy_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print('Exact adiabatic cooling for strain eps=1e-5:',exact,'K')
    for item in records:
        print(item['element'],item['thermal_scheme'],item['steps'],
            f"minT={item['min_temperature_K']:.6f}K",
            f"q_correction={item['mixed_timestep_coupling_correction_J']:.6f}J",
            f"relative_closure={item['relative_energy_identity_residual']:.2g}")
    return records,sign


if __name__=='__main__':
    analyze(Path(__file__).resolve().parents[2])
