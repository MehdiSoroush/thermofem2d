"""Reproduce curved circular-hole heat and isoparametric consistency benchmarks.

Outputs: CSV summary and a temperature-error convergence plot. This script
uses the verified circular geometry adapter; it does not change the ordinary
CLI geometry options or claim validated circular-hole stress concentrations.
"""
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from thermofem2d.circular_geometry import make_circular_mesh, integrated_geometry_area, hole_boundary_length
from thermofem2d.curved_field import circular_field_mesh
from thermofem2d.fem import assemble_thermal, assemble_elastic, boundary_load, solve_dirichlet, stresses

HERE = Path(__file__).resolve().parents[2]
CENTRE = np.array((0.12,-0.07),dtype=float)
RADIUS = 0.25


def measure(mat,kind,n):
    geo=make_circular_mesh(width=2.,height=2.,hole_radius=RADIUS,
                           hole_cx=float(CENTRE[0]),hole_cy=float(CENTRE[1]),
                           per_side=n,radial=n,element_type=kind)
    mesh=circular_field_mesh(geo)
    boundary=np.unique(np.concatenate(tuple(mesh.node_sets.values())))
    rad=np.linalg.norm(mesh.nodes-CENTRE,axis=1)
    exact_T=300+8*np.log(rad/RADIUS)
    K,C=assemble_thermal(mesh,mat)
    computed_T=solve_dirichlet(K,np.zeros(len(mesh.nodes)),
        {int(i):float(exact_T[i]) for i in boundary})
    rms_T=np.sqrt(np.mean((computed_T-exact_T)**2))
    # Constant strain and free thermal expansion are independent patch tests.
    A=np.array([[2.3e-5,-1.4e-5],[0.7e-5,1.1e-5]])
    imposed=(mesh.nodes@A.T+np.array([0.003,-0.001])).ravel()
    Tref=np.full(len(mesh.nodes),mat['T_ref'])
    Ku,f=assemble_elastic(mesh,mat,Tref)
    computed_u=solve_dirichlet(Ku,f,{2*int(i)+d:float(imposed[2*i+d]) for i in boundary for d in (0,1)})
    rms_u=np.sqrt(np.mean((computed_u-imposed)**2))
    return dict(element=kind,subdivisions=n,nodes=len(mesh.nodes),elements=len(mesh.elements),
       log_temperature_rms_error_K=float(rms_T),affine_displacement_rms_error_m=float(rms_u),
       integrated_area_m2=integrated_geometry_area(geo),
       plate_area_error_m2=abs(integrated_geometry_area(geo)-(4-np.pi*RADIUS**2)),
       circle_perimeter_error_m=abs(boundary_load(mesh,{'hole':1.0}).sum()-2*np.pi*RADIUS),
       thermal_capacity_constant_field_J_per_K=float(np.ones(len(mesh.nodes))@C@np.ones(len(mesh.nodes))))


def evaluate(output=HERE/'results/verification'):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    mat=json.loads((HERE/'examples/steady_plate.json').read_text())['material']
    records=[measure(mat,kind,n) for kind in ('T3','T6') for n in (3,4,6,8,10,12)]
    with (output/'circular_curved_field_convergence.csv').open('w',newline='',encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=list(records[0]),lineterminator='\n');writer.writeheader();writer.writerows(records)
    fig,ax=plt.subplots(figsize=(8,4.7))
    for kind in ('T3','T6'):
        dat=[x for x in records if x['element']==kind]
        ax.loglog([x['subdivisions'] for x in dat],
                  [x['log_temperature_rms_error_K'] for x in dat],'-o',label=kind)
    ax.set(xlabel='Subdivisions per plate side and radial direction',
           ylabel='RMS nodal temperature error (K)',
           title='Analytical harmonic temperature: eccentric circular hole')
    ax.grid(alpha=.25,which='both');ax.legend();fig.tight_layout()
    fig.savefig(output/'circular_curved_heat_convergence.png',dpi=175)
    plt.close(fig)
    print(f'Created {len(records)} physically computed benchmark rows.')
    for rec in records:
        print(f"{rec['element']} n={rec['subdivisions']:2d}: "
              f"temperature error={rec['log_temperature_rms_error_K']:.8g} K, "
              f"affine displacement error={rec['affine_displacement_rms_error_m']:.3g} m")
    return records


if __name__=='__main__':
    evaluate()
