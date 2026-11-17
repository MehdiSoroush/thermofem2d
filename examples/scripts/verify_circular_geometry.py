"""Reproduce analytic circular-hole geometry tests and mesh illustrations."""

from pathlib import Path
import csv
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from thermofem2d.circular_geometry import (
    make_circular_mesh, integrated_geometry_area, hole_boundary_length,
)


W = H = 2.0
R = 0.25
CX, CY = 0.12, -0.07


def make(kind, subdivisions):
    return make_circular_mesh(width=W, height=H, hole_radius=R,
                              hole_cx=CX, hole_cy=CY, per_side=subdivisions,
                              radial=subdivisions, element_type=kind)


def report(out=Path("results/verification")):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    expected_area = W*H-np.pi*R**2
    expected_perimeter = 2*np.pi*R
    for n in (3,4,6,8,12):
        for kind in ("T3", "T6"):
            geometry = make(kind, n)
            area = integrated_geometry_area(geometry)
            perimeter = hole_boundary_length(geometry)
            rows.append({"element":kind,"per_side":n,"radial":n,
                         "nodes":len(geometry.coordinates),
                         "triangles":len(geometry.triangles),
                         "area_m2":area,"area_error_m2":abs(area-expected_area),
                         "hole_length_m":perimeter,
                         "hole_length_error_m":abs(perimeter-expected_perimeter)})
    with (out/"circular_geometry_convergence.csv").open("w",newline="",encoding="utf-8") as fp:
        writer=csv.DictWriter(fp,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)

    fig,axs=plt.subplots(1,2,figsize=(11,4.5))
    for k,(field,ylabel) in enumerate((("area_error_m2","Absolute area error (m²)"),
                                      ("hole_length_error_m","Absolute hole-boundary length error (m)"))):
        for kind in ("T3","T6"):
            filtered=[r for r in rows if r["element"]==kind]
            axs[k].semilogy([r["per_side"] for r in filtered],
                            [r[field] for r in filtered],"o-",label=kind)
        axs[k].set(xlabel="Subdivisions per side",ylabel=ylabel)
        axs[k].grid(alpha=.25)
        axs[k].legend()
    fig.tight_layout();fig.savefig(out/"circular_geometry_convergence.png",dpi=180)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10,5))
    for ax,kind in zip(axes,("T3","T6")):
        g=make(kind,4)
        for tri in g.triangles:
            p=g.coordinates[tri[:3]]
            q=np.vstack((p,p[0]))
            ax.plot(q[:,0],q[:,1],linewidth=.5,color="0.5")
        for edge in g.boundary_edges["hole"]:
            p=g.coordinates[edge]
            if kind=="T3":
                xy=p
            else:
                xi=np.linspace(-1,1,15)
                N=np.column_stack((xi*(xi-1)/2,1-xi*xi,xi*(xi+1)/2))
                xy=N@p
            ax.plot(xy[:,0],xy[:,1],linewidth=1.2)
        phi=np.linspace(0,2*np.pi,400)
        ax.plot(CX+R*np.cos(phi),CY+R*np.sin(phi),"--",linewidth=.8,label="Exact circular boundary")
        ax.set(aspect="equal",title=f"{kind} · {len(g.triangles)} triangles",
               xlabel="x (m)",ylabel="y (m)")
        ax.legend(loc="lower right",fontsize=8)
    fig.tight_layout();fig.savefig(out/"circular_geometry_mesh.png",dpi=180)
    plt.close(fig)
    return rows


if __name__=="__main__":
    results=report()
    for item in results:
        print(item)
