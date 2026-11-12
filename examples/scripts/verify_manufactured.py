"""Reproduce exact-solution benchmarks for the T3 heat solver.

Run from repo root: python examples/scripts/verify_manufactured.py
No historical coursework results are used as numerical ground truth here.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from thermofem2d.fem import assemble_thermal, solve_dirichlet, element_heat_flux
from thermofem2d.geometry import make_mesh
from thermofem2d.solve import _heat_load


def benchmark(outdir=Path('results/verification')):
    root=Path(__file__).resolve().parents[2]
    config=json.loads((root/'examples/steady_plate.json').read_text(encoding='utf-8'))
    mat=config['material']
    out=Path(outdir);out.mkdir(parents=True, exist_ok=True)
    a=12.0
    rows=[]
    for n in (3, 4, 6, 8, 12, 16):
        geom={**config['geometry'],'per_side':n,'radial':n}
        mesh=make_mesh(**geom)
        xy=mesh.nodes
        exact=320+a*(xy[:,0]**2+xy[:,1]**2)
        boundary=np.unique(np.concatenate(tuple(mesh.node_sets.values())))
        fixed={int(i):float(exact[i]) for i in boundary}
        K,_=assemble_thermal(mesh,mat)
        f=_heat_load(mesh,mat,{'volumetric_heat_source':-4*a*mat['conductivity']})
        approx=solve_dirichlet(K,f,fixed)
        free=np.setdiff1d(np.arange(len(xy)),boundary)
        rms=float(np.linalg.norm((approx-exact)[free])/np.sqrt(len(free)))
        maxerr=float(np.max(np.abs(approx-exact)))
        rows.append({'subdivisions':n,'nodes':len(xy),'elements':len(mesh.elements),
                     'rms_free_node_error_K':rms,'max_nodal_error_K':maxerr})
    with (out/'quadratic_heat_convergence.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    ax.loglog([r['nodes'] for r in rows],[r['rms_free_node_error_K'] for r in rows],
              '-o',lw=2,label='T3 RMS free-node error')
    ax.set(xlabel='Number of mesh nodes',ylabel='Temperature RMS error [K]',
           title='Manufactured heat-conduction benchmark: T = 320 + 12(x² + y²)')
    ax.grid(True,which='both',alpha=.25);ax.legend()
    fig.savefig(out/'quadratic_heat_convergence.png',dpi=180);plt.close(fig)
    result={'verification':'analytical manufactured solution, not archival Abaqus',
            'exact_temperature':'T(x,y)=320+12*(x^2+y^2) K',
            'volumetric_source_W_per_m3':-4*a*mat['conductivity'],
            'temperature_dirichlet_boundary':'all outer and hole boundaries',
            'errors':rows,
            'observed_rms_error_reduction':rows[0]['rms_free_node_error_K']/rows[-1]['rms_free_node_error_K']}
    (out/'benchmark_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    benchmark()
