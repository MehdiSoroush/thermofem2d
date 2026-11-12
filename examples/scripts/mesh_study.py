"""Sensitivity to mesh density (NOT a singular-stress convergence proof).

Run: python examples/scripts/mesh_study.py
"""
import copy, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from thermofem2d.solve import solve_model

root=Path(__file__).resolve().parents[2]
cfg=json.loads((root/"examples/steady_plate.json").read_text())
rows=[]
for side,radial in [(4,4),(8,7),(12,10),(16,14),(20,17)]:
    case=copy.deepcopy(cfg)
    case["geometry"].update(per_side=side,radial=radial)
    r=solve_model(case)
    rows.append((len(r.mesh.nodes),len(r.mesh.elements),
                 float(np.linalg.norm(r.displacement.reshape(-1,2),axis=1).max()*1e3),
                 float(r.von_mises.max()/1e6)))
out=root/"results"; (out/"data").mkdir(parents=True,exist_ok=True)
with (out/"data"/"mesh_sensitivity.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["nodes","elements","max_displacement_mm","max_centroid_von_mises_MPa"]);w.writerows(rows)
fig,ax=plt.subplots(figsize=(8,5),layout="constrained")
ax.plot([x[1] for x in rows],[x[2] for x in rows],"o-",label="Maximum displacement")
ax.set(xlabel="T3 element count",ylabel="Maximum displacement [mm]",title="Mesh sensitivity | displacement")
ax.grid(alpha=.25)
p=out/"figures"/"mesh_sensitivity_displacement.png";p.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(p,dpi=180);plt.close(fig)
print("Mesh-sensitivity results:")
for x in rows:print(x)
print("NOTE: peak centroid stress continues to change; do not claim stress convergence.")
