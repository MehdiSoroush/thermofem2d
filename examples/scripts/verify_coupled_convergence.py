"""v0.4 reproducible temporal/energy and exploratory spatial sensitivity checks.

Reference: independent DOP853 integration of the *same* semi-discrete FEM ODE.
Results are not external FE validation or proof of convergence of peak stress.
"""
import copy
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from thermofem2d.solve import solve_model
from thermofem2d.verification import solve_semidiscrete_ode,insulated_energy_audit


def case(kind, subdiv, steps):
    cfg=json.loads((ROOT/'examples/coupled_plate_t6.json').read_text(encoding='utf8'))
    cfg['geometry'].update(element_type=kind,per_side=subdiv,radial=subdiv)
    cfg['analysis'].update(steps=int(steps),time_step=1e-4/steps,
        initial_temperature={'type':'affine','base':310.,'gx':6.,'gy':-3.},
        temperature_bc={},traction_bc={},inward_heat_flux={},
        volumetric_heat_source=0.,thermoelastic_feedback=True)
    return cfg


def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def run(out=None):
    out=Path(out) if out is not None else ROOT/'results'/'verification'
    out.mkdir(parents=True,exist_ok=True)
    rows=[];energies=[]
    for kind in ('T3','T6'):
        # Identical spatial mesh for dt refinement; no spatial error cancellation claim.
        for steps in (8,16,32,64):
            cfg=case(kind,2,steps)
            sol=solve_model(cfg)
            ref=solve_semidiscrete_ode(cfg,sol.times)
            en=insulated_energy_audit(cfg,sol)
            initial_energy=en.mechanical_energy[0]+en.thermal_energy[0]
            rows.append(dict(element=kind,elements=len(sol.mesh.elements),nodes=len(sol.mesh.nodes),
                dt_s=cfg['analysis']['time_step'],steps=steps,
                max_final_temperature_error_K=float(np.max(abs(sol.temperature-ref.temperature_history[-1]))),
                max_final_displacement_error_m=float(np.max(abs(sol.displacement-ref.displacement_history[-1]))),
                max_final_velocity_error_m_s=float(np.max(abs(sol.velocity_history[-1]-ref.velocity_history[-1]))),
                ode_residual=ref.max_relative_residual))
            energies.append(dict(element=kind,dt_s=cfg['analysis']['time_step'],steps=steps,
                initial_energy_J=initial_energy,
                final_energy_J=float(en.mechanical_energy[-1]+en.thermal_energy[-1]),
                conduction_loss_J=float(en.cumulative_conduction[-1]),
                final_relative_energy_balance_error=float(en.energy_balance_error[-1]/initial_energy),
                max_absolute_relative_energy_balance_error=float(np.max(abs(en.energy_balance_error))/initial_energy)))
    write_csv(out/'coupled_time_convergence.csv',rows)
    write_csv(out/'coupled_energy_balance.csv',energies)

    fig,ax=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for kind in ('T3','T6'):
        r=[x for x in rows if x['element']==kind]
        x=np.array([i['dt_s'] for i in r]);
        ax[0].loglog(x,[i['max_final_temperature_error_K'] for i in r],'-o',label=kind)
        ax[1].loglog(x,[i['max_final_displacement_error_m']*1e6 for i in r],'-o',label=kind)
    for axis in ax:
        axis.set_xlabel('Time step [s]');axis.grid(alpha=.2,which='both');axis.legend()
    ax[0].set_ylabel('Max final temperature difference [K]')
    ax[1].set_ylabel('Max final displacement difference [µm]')
    fig.suptitle('Monolithic Newmark/BE vs independent DOP853 time integration\n(same semi-discrete T3/T6 matrices)')
    fig.savefig(out/'coupled_time_convergence.png',dpi=180);plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.2,4.3),layout='constrained')
    for kind in ('T3','T6'):
        r=[x for x in energies if x['element']==kind]
        ax.loglog([i['dt_s'] for i in r],[abs(i['final_relative_energy_balance_error']) for i in r],'-o',label=kind)
    ax.set(xlabel='Time step [s]',ylabel='|Final energy-balance error| / E₀',
           title='Insulated source-free coupled system: energy-budget defect')
    ax.grid(alpha=.2,which='both');ax.legend()
    fig.savefig(out/'coupled_energy_balance.png',dpi=180);plt.close(fig)

    spatial=[]
    for kind in ('T3','T6'):
        for subdiv in (2,3,4,6):
            cfg=case(kind,subdiv,32)
            tic=time.perf_counter();sol=solve_model(cfg);elapsed=time.perf_counter()-tic
            # Mean temperature is area weighted, not a raw nodal arithmetic mean.
            vertices=sol.mesh.nodes[sol.mesh.elements[:,:3]]
            area=.5*np.abs((vertices[:,1,0]-vertices[:,0,0])*(vertices[:,2,1]-vertices[:,0,1])-(vertices[:,1,1]-vertices[:,0,1])*(vertices[:,2,0]-vertices[:,0,0]))
            conn=sol.mesh.elements
            if kind=='T3': mean_field=sol.temperature[conn].mean(axis=1)
            else: mean_field=sol.temperature[conn[:,3:]].mean(axis=1) # exact degree-2 integration, three edge midpoints
            avg=float(area@mean_field/area.sum())
            spatial.append(dict(element=kind,elements=len(conn),nodes=len(sol.mesh.nodes),
                mechanical_dofs=2*len(sol.mesh.nodes),thermal_dofs=len(sol.mesh.nodes),
                area_avg_temperature_K=avg,max_displacement_m=float(np.linalg.norm(sol.displacement.reshape(-1,2),axis=1).max()),
                max_centroid_von_mises_Pa=float(sol.von_mises.max()),elapsed_wall_s=elapsed))
    write_csv(out/'coupled_spatial_sensitivity.csv',spatial)
    fig,ax=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for kind in ('T3','T6'):
        r=[x for x in spatial if x['element']==kind]
        ax[0].plot([i['mechanical_dofs'] for i in r],[i['max_displacement_m']*1e6 for i in r],'-o',label=kind)
        ax[1].plot([i['mechanical_dofs'] for i in r],[i['max_centroid_von_mises_Pa']/1e6 for i in r],'-o',label=kind)
    for axis in ax: axis.set(xlabel='Mechanical degrees of freedom');axis.legend();axis.grid(alpha=.2)
    ax[0].set_ylabel('Maximum nodal displacement [µm]');ax[1].set_ylabel('Maximum centroid von Mises [MPa]')
    fig.suptitle('Exploratory mesh sensitivity; no continuum exact solution asserted')
    fig.savefig(out/'coupled_spatial_sensitivity.png',dpi=180);plt.close(fig)
    print(json.dumps({'time_convergence_rows':len(rows),'energy_rows':len(energies),'spatial_rows':len(spatial),
        'max_ode_residual':max(x['ode_residual'] for x in rows),
        'largest_relative_energy_error':max(abs(x['final_relative_energy_balance_error']) for x in energies)},indent=2))
    return rows,energies,spatial

if __name__=='__main__':
    run()
