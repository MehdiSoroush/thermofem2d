"""Optional interactive explorer. Requires `pip install -e '.[app]'`."""
from __future__ import annotations
import copy,io,json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import streamlit as st
from thermofem2d.solve import solve_model
from thermofem2d.element import plotting_triangles,element_type_for_mesh

st.set_page_config(page_title="ThermoFEM-2D | FEM laboratory", layout="wide")
st.title("ThermoFEM-2D")
st.caption("v0.5.0 • Parametric T3/T6 plate FEM • sequential or fully coupled transient dynamics")
st.warning("Research prototype. The coupled option is a new reduced plane-stress T3/T6 implementation; it has not been benchmarked against the original 2016 numerical results or independent industrial software.")
base=json.loads((Path(__file__).parent/"examples/steady_plate.json").read_text())
config=copy.deepcopy(base)
with st.sidebar:
    st.header("Geometry [m]")
    g=config["geometry"]
    g["hole_width"]=st.slider("Hole width",0.10,0.90,float(g["hole_width"]),0.02)
    g["hole_height"]=st.slider("Hole height",0.10,0.90,float(g["hole_height"]),0.02)
    g["hole_cx"]=st.slider("Hole x-position",-0.45,0.45,0.0,0.02)
    g["hole_cy"]=st.slider("Hole y-position",-0.45,0.45,0.0,0.02)
    g["per_side"]=st.slider("Boundary subdivisions per side",4,28,12)
    g["radial"]=st.slider("Radial subdivisions",4,24,12)
    g["element_type"]=st.selectbox("Finite element",["T3","T6"],index=0)
    st.header("Thermal loading [K]")
    a=config["analysis"]
    a["temperature_bc"]["left"]=st.number_input("Left boundary temperature",value=380.0)
    a["temperature_bc"]["right"]=st.number_input("Right boundary temperature",value=298.0)
    st.header("Physics")
    a["thermoelastic_mode"] = st.selectbox("Thermomechanical physics",
        ["one_way_quasistatic", "coupled_dynamic"],
        format_func=lambda x: "Sequential quasi-static (v0.2 baseline)" if x=="one_way_quasistatic" else "Two-way coupled dynamic (v0.3)")
    if a["thermoelastic_mode"] == "coupled_dynamic":
        a["thermal_mode"] = "transient"
        a["thermoelastic_feedback"] = st.checkbox("Enable strain-rate → temperature feedback", value=True)
        a["initial_temperature"] = st.number_input("Initial interior temperature [K]", value=298.0)
        a["time_step"] = st.number_input("Time step [s]", min_value=1e-8, value=4e-5, format="%.8f")
        a["steps"] = st.slider("Dynamic time steps", 1, 60, 12)
        st.caption("Constant temperature and displacement boundaries; idealized small-strain plane stress.")
    else:
        a["thermal_mode"]=st.selectbox("Thermal model",["steady","transient"])
        if a["thermal_mode"]=="transient":
            a["initial_temperature"]=st.number_input("Initial interior temperature [K]",value=298.0)
            a["time_step"]=st.number_input("Time step [s]",min_value=1.0,value=14400.0)
            a["steps"]=st.slider("Time steps",1,120,10)
    run=st.button("Solve model",type="primary",use_container_width=True)

st.markdown("**Model assumptions:** homogeneous isotropic material, constant properties, Fourier conduction, small-strain plane stress, traction-free hole by default. Two-way coupling is optional and uses a linearized constant-reference-temperature model.")
if not run:
    st.info("Configure the model in the sidebar and select **Solve model**.")
else:
    try:
        with st.spinner("Assembling and solving FEM equations"):
            sol=solve_model(config)
    except (ValueError,KeyError,RuntimeError) as exc:
        st.error(f"Model could not be solved: {exc}")
    else:
        xy=sol.mesh.nodes
        kind=element_type_for_mesh(sol.mesh)
        subdivisions=4 if kind=="T6" else 1
        tri=mtri.Triangulation(xy[:,0],xy[:,1],plotting_triangles(sol.mesh))
        st.success(f"Solved: {len(xy)} nodes, {len(sol.mesh.elements)} {kind} elements")
        cols=st.columns(3)
        cols[0].metric("Maximum temperature",f"{sol.temperature.max():.1f} K")
        cols[1].metric("Maximum displacement",f"{np.linalg.norm(sol.displacement.reshape(-1,2),axis=1).max()*1000:.4g} mm")
        cols[2].metric("Maximum centroid von Mises",f"{sol.von_mises.max()/1e6:.3g} MPa")
        cases=[("Temperature [K]",sol.temperature,"inferno",False),
               ("Displacement [mm]",np.linalg.norm(sol.displacement.reshape(-1,2),axis=1)*1e3,"viridis",False),
               ("Element-centroid von Mises [MPa]",sol.von_mises/1e6,"magma",True)]
        for col,(title,field,cmap,elemental) in zip(st.columns(3),cases):
            fig,ax=plt.subplots(figsize=(5.5,4.4),layout="constrained")
            if elemental: im=ax.tripcolor(tri,facecolors=np.repeat(field,subdivisions),shading="flat",cmap=cmap)
            else: im=ax.tripcolor(tri,field,shading="gouraud",cmap=cmap)
            ax.set_aspect("equal");ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.75)
            col.pyplot(fig,clear_figure=True)
            plt.close(fig)
        with st.expander("Inspect mesh and time response"):
            fig,ax=plt.subplots(figsize=(6,5))
            ax.triplot(tri,lw=.4,color="#19857D");ax.set_aspect("equal")
            st.pyplot(fig);plt.close(fig)
            if len(sol.times)>1:
                if sol.displacement_history is not None:
                    st.line_chart({"Time [s]":sol.times,
                                   "Mean nodal temperature [K]":sol.temperature_history.mean(axis=1)},
                                  x="Time [s]", y="Mean nodal temperature [K]")
                    maximum_u = np.linalg.norm(sol.displacement_history.reshape(len(sol.times),-1,2),axis=2).max(axis=1)
                    st.line_chart({"Time [s]":sol.times, "Max displacement [mm]":1000*maximum_u},
                                  x="Time [s]",y="Max displacement [mm]")
                else:
                    st.line_chart({"Time [h]":sol.times/3600,
                                   "Mean nodal temperature [K]":sol.temperature_history.mean(axis=1)},x="Time [h]",y="Mean nodal temperature [K]")
        st.download_button("Download model JSON",data=json.dumps(config,indent=2),file_name="thermofem2d_model.json",mime="application/json")
        buf=io.StringIO();buf.write("x_m,y_m,T_K,ux_m,uy_m\n")
        np.savetxt(buf,np.column_stack((xy,sol.temperature,sol.displacement[::2],sol.displacement[1::2])),delimiter=",")
        st.download_button("Download nodal data CSV",data=buf.getvalue(),file_name="nodal_fields.csv",mime="text/csv")
