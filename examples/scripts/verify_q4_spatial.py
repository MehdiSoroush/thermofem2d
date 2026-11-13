"""Reproduce independent Q4/T3/T6 spatial comparisons and affine patch checks."""
from pathlib import Path
import copy
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

from thermofem2d.q4_reference import quad_mesh,q4_patch_tests,solve_q4_coupled
from thermofem2d.solve import solve_model


def evaluate(out=Path('results/verification')):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    base=json.loads(Path('examples/coupled_plate_t6.json').read_text(encoding='utf-8'))
    results=[];patch=[]
    for n in (3,4,6,8):
        c=copy.deepcopy(base)
        c['geometry'].update(per_side=n,radial=n)
        c['analysis'].update(steps=8,time_step=2e-5)
        q4=solve_q4_coupled(c)
        terr, uerr=q4_patch_tests(quad_mesh(c['geometry']),c['material'])
        patch.append(dict(per_side=n,q4_elements=len(q4.mesh.cells),temperature_patch_error_K=terr,
                          displacement_patch_error_m=uerr))
        for kind in ('T3','T6'):
            c['geometry']['element_type']=kind
            tri=solve_model(c)
            dist,idx=cKDTree(tri.mesh.nodes).query(q4.mesh.nodes)
            if dist.max()>1e-12:
                raise ValueError('Q4 reference nodes do not match production mesh vertices')
            Tdiff=tri.temperature[idx]-q4.temperature
            Udiff=tri.displacement.reshape(-1,2)[idx]-q4.displacement.reshape(-1,2)
            errT=float(np.sqrt(np.mean(Tdiff*Tdiff)))
            errU=float(np.sqrt(np.mean(np.sum(Udiff*Udiff,axis=1))))
            ref_mag=float(np.sqrt(np.mean(np.sum(q4.displacement.reshape(-1,2)**2,axis=1))))
            results.append(dict(element=kind,per_side=n,parent_triangles=len(tri.mesh.elements),
                                q4_elements=len(q4.mesh.cells),shared_vertices=len(q4.mesh.nodes),
                                rms_temperature_difference_K=errT,
                                rms_displacement_vector_difference_m=errU,
                                relative_displacement_difference=errU/max(ref_mag,1e-18)))
        print('Compared refinement',n,'Q4 patch errors',terr,uerr)
    for name,rows in [('q4_spatial_comparison.csv',results),('q4_patch_tests.csv',patch)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for key,ylabel,filename in [('rms_temperature_difference_K','RMS final-temperature difference (K)','q4_temperature_comparison.png'),
                                 ('rms_displacement_vector_difference_m','RMS final-displacement difference (m)','q4_displacement_comparison.png')]:
        fig,ax=plt.subplots(figsize=(7,4.4))
        for kind in ('T3','T6'):
            r=[item for item in results if item['element']==kind]
            ax.plot([z['per_side'] for z in r],[z[key] for z in r],'-o',label=f'{kind} vs independent Q4')
        ax.set(xlabel='Subdivisions per side = radial layers',ylabel=ylabel)
        ax.set_yscale('log');ax.grid(True,alpha=.25);ax.legend();fig.tight_layout()
        fig.savefig(out/filename,dpi=170);plt.close(fig)
    return results,patch


if __name__=='__main__':
    evaluate()
