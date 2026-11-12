"""Command-line entry point; creates auditable figures and data outputs."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from .solve import solve_model
from .plot import make_figures
from .fem import element_heat_flux


def main(argv=None):
    p=argparse.ArgumentParser(description="ThermoFEM-2D: perforated-plate FEM prototype")
    p.add_argument("--config",default="examples/steady_plate.json",help="JSON model input")
    p.add_argument("--out",default="results",help="output directory")
    args=p.parse_args(argv)
    config=json.loads(Path(args.config).read_text(encoding="utf-8"))
    sol=solve_model(config)
    out=Path(args.out);(out/"data").mkdir(parents=True,exist_ok=True)
    images=make_figures(sol,out/"figures")
    n=len(sol.mesh.nodes)
    np.savetxt(out/"data"/"nodal_fields.csv",
               np.column_stack((sol.mesh.nodes,sol.temperature,sol.displacement[::2],sol.displacement[1::2])),
               delimiter=",",header="x_m,y_m,temperature_K,ux_m,uy_m",comments="")
    np.savetxt(out/"data"/"element_stresses.csv",
               np.column_stack((sol.sigma,sol.von_mises)),
               delimiter=",",header="sigma_xx_Pa,sigma_yy_Pa,tau_xy_Pa,von_mises_Pa",comments="")
    flux = element_heat_flux(sol.mesh, config["material"], sol.temperature)
    centroids = sol.mesh.nodes[sol.mesh.elements].mean(axis=1)
    np.savetxt(out/"data"/"element_heat_flux.csv",
               np.column_stack((centroids, flux)),
               delimiter=",",header="centroid_x_m,centroid_y_m,heat_flux_x_W_per_m2,heat_flux_y_W_per_m2",comments="")
    nodes=sol.mesh.nodes
    tri=sol.mesh.elements
    p0,p1,p2=nodes[tri[:,0]],nodes[tri[:,1]],nodes[tri[:,2]]
    material_area = 0.5 * np.sum((p1[:,0]-p0[:,0])*(p2[:,1]-p0[:,1]) - (p1[:,1]-p0[:,1])*(p2[:,0]-p0[:,0]))
    summary={"version":"0.1.1","source_config":str(args.config),
             "nodes":n,
             "plate_material_area_m2":float(material_area),
             "max_element_heat_flux_W_per_m2":float(np.linalg.norm(flux,axis=1).max()),"elements":len(sol.mesh.elements),
             "temperature_min_K":float(sol.temperature.min()),
             "temperature_max_K":float(sol.temperature.max()),
             "max_displacement_m":float(np.max(np.linalg.norm(sol.displacement.reshape(-1,2),axis=1))),
             "max_element_centroid_von_mises_Pa":float(sol.von_mises.max()),
             "analysis_type":"one-way thermal to quasi-static plane stress",
             "validation_status":"analytical and manufactured-solution tests; historical coupled model not reproduced"}
    (out/"data"/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2));print("Figures:",*[str(x) for x in images],sep="\n  ")
    return summary


if __name__=="__main__":
    main()
