"""Reproducible plots; stress and flux are reported at parent-element centroids."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from .solve import Solution
from .fem import element_heat_flux
from .element import plotting_triangles,element_type_for_mesh


def make_figures(s:Solution,outdir:str|Path)->list[Path]:
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    xy=s.mesh.nodes; parents=s.mesh.elements
    element_kind=element_type_for_mesh(s.mesh)
    triangles=plotting_triangles(s.mesh)
    refinement=4 if element_kind=="T6" else 1
    triang=mtri.Triangulation(xy[:,0],xy[:,1],triangles)
    N=len(xy)
    ux=s.displacement[0:2*N:2]
    uy=s.displacement[1:2*N:2]
    magnitude=np.hypot(ux,uy)*1e3
    fields=[("temperature",s.temperature,"Temperature [K]","inferno","node"),
            ("displacement",magnitude,"Displacement magnitude [mm]","viridis","node"),
            ("von_mises",s.von_mises/1e6,"Von Mises stress [MPa]","magma","element")]
    output=[]
    for name,data,label,cmap,field_kind in fields:
        fig,ax=plt.subplots(figsize=(8,6),layout="constrained")
        if field_kind=="element":
            plot=ax.tripcolor(triang,facecolors=np.repeat(data,refinement),shading="flat",cmap=cmap)
        else:
            plot=ax.tripcolor(triang,data,shading="gouraud",cmap=cmap)
        ax.set_aspect("equal");ax.set_xlabel("x [m]");ax.set_ylabel("y [m]")
        ax.set_title(label+f" | ThermoFEM-2D {element_kind}")
        fig.colorbar(plot,ax=ax,label=label,shrink=.86)
        p=outdir/(name+".png");fig.savefig(p,dpi=180);plt.close(fig);output.append(p)
    fig,ax=plt.subplots(figsize=(8,6),layout="constrained")
    ax.triplot(triang,lw=.46,color="#137c84",alpha=.65)
    ax.scatter(xy[:,0],xy[:,1],s=1,c="#24364a")
    ax.set_aspect("equal");ax.set_xlabel("x [m]");ax.set_ylabel("y [m]")
    ax.set_title(f"Geometry-fitted triangular mesh | {len(xy)} nodes, {len(parents)} {element_kind} elements")
    p=outdir/"mesh.png";fig.savefig(p,dpi=180);plt.close(fig);output.append(p)
    # Raw elementwise Fourier flux; color is magnitude, arrows are direction.
    flux=element_heat_flux(s.mesh, s.material, s.temperature)
    centers=xy[parents[:,:3]].mean(axis=1)
    fluxmag=np.linalg.norm(flux,axis=1)
    fig,ax=plt.subplots(figsize=(8,6),layout="constrained")
    im=ax.tripcolor(triang,facecolors=np.repeat(fluxmag,refinement),shading="flat",cmap="cividis")
    picks=np.arange(0,len(flux),max(1,int(np.ceil(len(flux)/120))))
    qscale=max(float(fluxmag.max()),1e-30)
    norm_len=0.045*max(float(np.ptp(xy[:,0])),float(np.ptp(xy[:,1])))
    vec=(flux[picks]/qscale)*norm_len
    ax.quiver(centers[picks,0],centers[picks,1],vec[:,0],vec[:,1],
              color="white",angles="xy",scale_units="xy",scale=1,
              width=0.003,headwidth=3.4,headlength=4.5)
    ax.set(aspect="equal",xlabel="x [m]",ylabel="y [m]",
           title=f"Fourier heat flux magnitude and direction | {element_kind} centroids")
    fig.colorbar(im,ax=ax,label="Heat flux magnitude [W/m²]",shrink=.86)
    p=outdir/"heat_flux.png";fig.savefig(p,dpi=180);plt.close(fig);output.append(p)

    # Original and magnified deformed meshes. Always state scale factor.
    u=s.displacement.reshape(-1,2)
    umax=float(np.linalg.norm(u,axis=1).max())
    length=float(max(np.ptp(xy[:,0]),np.ptp(xy[:,1])))
    scale=0.08*length/umax if umax>0 else 1.0
    xyz=xy+scale*u
    deformation=mtri.Triangulation(xyz[:,0],xyz[:,1],triangles)
    fig,ax=plt.subplots(figsize=(8,6),layout="constrained")
    ax.triplot(triang,color="gray",lw=.25,alpha=.3,label="Undeformed")
    ax.triplot(deformation,color="#087f8c",lw=.45,alpha=.70,label="Deformed (magnified)")
    ax.set(aspect="equal",xlabel="x [m]",ylabel="y [m]",
           title=f"Original vs deformed mesh | displacement ×{scale:.4g}")
    ax.legend(loc="upper right")
    p=outdir/"deformed_mesh.png";fig.savefig(p,dpi=180);plt.close(fig);output.append(p)
    fig,axes=plt.subplots(2,3,figsize=(16,10),layout="constrained")
    for ax,(name,data,label,cmap,field_kind) in zip(axes[0],fields):
        if field_kind=="element":
            im=ax.tripcolor(triang,facecolors=np.repeat(data,refinement),shading="flat",cmap=cmap)
        else:
            im=ax.tripcolor(triang,data,shading="gouraud",cmap=cmap)
        ax.set_aspect("equal");ax.set_xlabel("x [m]");ax.set_ylabel("y [m]")
        ax.set_title(label);fig.colorbar(im,ax=ax,shrink=.80)
    ax=axes[1,0]
    im=ax.tripcolor(triang,facecolors=np.repeat(fluxmag,refinement),shading="flat",cmap="cividis")
    ax.quiver(centers[picks,0],centers[picks,1],vec[:,0],vec[:,1],
              color="white",angles="xy",scale_units="xy",scale=1,width=.003)
    ax.set_title("Heat flux [W/m²] | elementwise");fig.colorbar(im,ax=ax,shrink=.80)
    ax=axes[1,1];ax.triplot(triang,lw=.38,color="#137c84",alpha=.65)
    ax.set_title(f"{element_kind} mesh | {len(parents)} parent elements")
    ax=axes[1,2]
    ax.triplot(triang,color="gray",lw=.25,alpha=.30)
    ax.triplot(deformation,color="#087f8c",lw=.40,alpha=.70)
    ax.set_title(f"Deformation displayed ×{scale:.3g}")
    for ax in axes[1]:
        ax.set_aspect("equal");ax.set_xlabel("x [m]");ax.set_ylabel("y [m]")
    fig.suptitle(f"ThermoFEM-2D v0.2.0 | {element_kind} computed fields and mesh")
    p=outdir/"overview.png";fig.savefig(p,dpi=180);plt.close(fig);output.append(p)
    if len(s.times)>1:
        fig,ax=plt.subplots(figsize=(9,5),layout="constrained")
        # Spatial mean of nodal temperatures (diagnostic, not mass-weighted).
        ax.plot(s.times/3600,np.mean(s.temperature_history,axis=1),lw=2)
        ax.set(xlabel="Time [h]",ylabel="Mean nodal temperature [K]",
               title="Transient thermal response | backward Euler")
        ax.grid(alpha=.2)
        p=outdir/"thermal_history.png";fig.savefig(p,dpi=180);plt.close(fig);output.append(p)
    return output
