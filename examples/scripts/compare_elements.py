"""Reproducible T3/T6 comparison on the perforated-plate demonstration case.

CAUTION: this is a sensitivity study, not an independent reference solution.
Element-centroid peak stresses are *not* guaranteed mesh convergent.
"""
import csv,copy,json,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from thermofem2d.solve import solve_model


def run_comparison(output:Path|None=None):
    cfg=json.loads((ROOT/'examples'/'steady_plate.json').read_text())
    out=output or ROOT/'results'/'verification'
    out.mkdir(parents=True,exist_ok=True)
    records=[]
    for subdivisions in (4,6,8,10):
        for element_type in ('T3','T6'):
            case=copy.deepcopy(cfg)
            case['geometry'].update(per_side=subdivisions,radial=subdivisions,
                                    element_type=element_type)
            s=solve_model(case)
            records.append(dict(element_type=element_type,
                                per_side=subdivisions,radial=subdivisions,
                                elements=len(s.mesh.elements),nodes=len(s.mesh.nodes),
                                thermal_dofs=len(s.mesh.nodes),
                                mechanical_dofs=2*len(s.mesh.nodes),
                                min_temperature_K=float(s.temperature.min()),
                                max_temperature_K=float(s.temperature.max()),
                                max_displacement_mm=float(np.linalg.norm(s.displacement.reshape(-1,2),axis=1).max()*1000),
                                max_centroid_von_mises_MPa=float(s.von_mises.max()/1e6)))
    csvfile=out/'t3_t6_comparison.csv'
    with csvfile.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=records[0]);writer.writeheader();writer.writerows(records)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    for element_type in ('T3','T6'):
        a=[r for r in records if r['element_type']==element_type]
        x=[r['mechanical_dofs'] for r in a]
        axes[0].plot(x,[r['max_displacement_mm'] for r in a],'-o',label=element_type)
        axes[1].plot(x,[r['max_centroid_von_mises_MPa'] for r in a],'-o',label=element_type)
    axes[0].set(ylabel='Max nodal displacement [mm]',xlabel='Mechanical DOFs',title='Global displacement sensitivity')
    axes[1].set(ylabel='Max element-centroid von Mises [MPa]',xlabel='Mechanical DOFs',title='Stress peaks not certified convergent')
    for ax in axes:
        ax.grid(alpha=.2);ax.legend()
    fig.suptitle('T3 versus affine T6 | one-way plane-stress demonstration')
    fig.savefig(out/'t3_t6_comparison.png',dpi=180)
    plt.close(fig)
    return records


if __name__=='__main__':
    for row in run_comparison():print(row)
