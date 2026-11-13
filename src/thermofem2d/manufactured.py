"""Continuum-manufactured, fully coupled transient thermoelastic benchmark.

The analytic displacement/temperature fields satisfy the *strong continuum PDE*
with independently evaluated body force and volumetric heat source. No solver
or FEM assembly routines are used to construct these forcing functions.

For the verification harness, all boundary nodes follow the exact time-varying
Dirichlet fields. The production public API does not yet support these BCs.
The harness uses the production T3/T6 or independently assembled Q4 matrices.
It tests spatial approximation, coupling blocks and the time stepper together.

This is NOT external-package validation, experimental validation, or a proof
that a sharp-corner singular stress is converged.
"""
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.sparse import bmat
from scipy.sparse.linalg import splu

from .element import TRIANGLE_QUADRATURE, shape_and_gradients
from .geometry import make_mesh
from .fem import assemble_elastic, assemble_thermal
from .coupled import assemble_inertia_coupling
from .q4_reference import quad_mesh, assemble_quad, q4_shape, _GAUSS


@dataclass
class ManufacturedResult:
    element: str
    subdivisions: int
    nodes: int
    elements: int
    steps: int
    temperature_rms_error_K: float
    displacement_rms_error_m: float
    temperature_max_error_K: float
    displacement_max_error_m: float
    thermal_scale_K: float
    displacement_scale_m: float
    solved_temperature: np.ndarray
    solved_displacement: np.ndarray
    exact_temperature: np.ndarray
    exact_displacement: np.ndarray


def exact_fields(xy: np.ndarray, time: float, *, t_ref=298.0, a=0.02, b=6000.0, curvature=0.0):
    """u=a*(t+curvature*t²)*(x²,xy), theta=b*t*(x²+2y²).

    Both fields are spatially quadratic. Nonzero curvature gives constant
    mechanical acceleration and linearly varying strain-rate feedback.
    """
    x, y = xy[:, 0], xy[:, 1]
    g=a*(time+curvature*time*time)
    return (np.column_stack((g*x*x,g*x*y)).ravel(),
            t_ref + b * time * (x*x + 2*y*y))


def exact_initial_velocity(xy: np.ndarray, a=0.02):
    x,y=xy.T
    return np.column_stack((a*x*x, a*x*y)).ravel()


def continuum_sources(xy: np.ndarray, time: float, material, *, a=0.02, b=6000.0, curvature=0.0):
    """Analytically derived (f_x,f_y) [N/m³], q [W/m³].

    Strong PDE with plane-stress reduction (theta=T-Tref):
      rho*u_tt - div(D*eps(u) - beta*theta*I) = f
      rho*c*theta_t - k*laplacian(theta) + Tref*beta*div(u_t) = q
    beta=E*alpha/(1-nu); no gravity, boundary flux, or other sources.

    In the reference fields, u_tt=2*a*curvature*(x²,xy),
    laplace(theta)=6*b*t, div(u_t)=3*a*(1+2*curvature*t)*x,
    and div(D eps(u))=[(C(2+nu)+mu)*a*(t+curvature*t²),0].
    """
    x,y=xy[:,0],xy[:,1]
    E=float(material['E']);nu=float(material['nu'])
    alpha=float(material['alpha']);rho=float(material['density'])
    cap=float(material['heat_capacity']);k=float(material['conductivity'])
    tref=float(material['T_ref'])
    C=E/(1-nu**2);mu=E/(2*(1+nu));beta=E*alpha/(1-nu)
    g=a*(time+curvature*time*time)
    inertial=2*a*curvature*rho
    fx=inertial*x*x-(C*(2+nu)+mu)*g + beta*(2*b*time*x)
    fy=inertial*x*y+beta*(4*b*time*y)
    q=rho*cap*b*(x*x+2*y*y)-k*(6*b*time)+tref*beta*(3*a*(1+2*curvature*time)*x)
    return np.column_stack((fx,fy)),q


def _loads(mesh, kind, material, time, curvature):
    """Integrate analytic body/heat sources against FEM N; not via FEM residual."""
    n=len(mesh.nodes);h=material['thickness']
    f=np.zeros(2*n);q=np.zeros(n)
    if kind=='Q4':
        gauss=[(xi,eta,1.0) for xi in _GAUSS for eta in _GAUSS]
        for conn in mesh.cells:
            coords=mesh.nodes[conn]
            for xi,eta,w in gauss:
                N,ref_grad=q4_shape(xi,eta)
                det=float(np.linalg.det(coords.T@ref_grad))
                if det<=0: raise ValueError('inverted Q4 mesh')
                gp=N@coords
                forces,heat=continuum_sources(gp[None,:],time,material,curvature=curvature)
                fac=h*w*det
                f[2*conn] +=fac*N*forces[0,0]
                f[2*conn+1]+=fac*N*forces[0,1]
                q[conn] +=fac*N*heat[0]
    else:
        for conn in mesh.elements:
            coords=mesh.nodes[conn[:3]]
            for r,s,w in TRIANGLE_QUADRATURE:
                N,_,area=shape_and_gradients(coords,kind,r,s)
                gp=N@mesh.nodes[conn]
                forces,heat=continuum_sources(gp[None,:],time,material,curvature=curvature)
                fac=h*w*area
                f[2*conn]+=fac*N*forces[0,0]
                f[2*conn+1]+=fac*N*forces[0,1]
                q[conn]+=fac*N*heat[0]
    return f,q


def verify_manufactured_coupled(config, kind:Literal['T3','T6','Q4'],
                                subdivisions:int=4, steps:int=8, end_time:float=1e-3,
                                *, disable_feedback:bool=False, curvature:float=0.0):
    """Solve independent continuum MMS with prescribed exact moving boundaries.

    Matrices use production T3/T6 or separately assembled Q4. Loads and exact
    fields are computed from analytic continuum formulae. The displacement
    and temperature are initialized exactly, including nonzero initial velocity.
    Fixed velocity/acceleration track the exact polynomial-time BCs.
    With disable_feedback=True, the source remains the *fully coupled* exact
    source, providing a deliberate negative control for this benchmark.
    """
    if kind not in ('T3','T6','Q4'):raise ValueError('kind must be T3, T6 or Q4')
    if not np.isfinite(curvature):raise ValueError('curvature must be finite')
    if isinstance(steps,bool) or not isinstance(steps,int) or steps<1:raise ValueError('steps >=1')
    if isinstance(subdivisions,bool) or not isinstance(subdivisions,int) or subdivisions<3:raise ValueError('subdivisions >=3')
    dt=end_time/steps
    if not np.isfinite(dt) or dt<=0:raise ValueError('end_time must be finite positive')
    g=dict(config['geometry'],per_side=subdivisions,radial=subdivisions)
    mat=config['material'];tref=mat['T_ref']
    if kind=='Q4':
        mesh=quad_mesh(g)
        Kt,Ct,Ku,M,G=assemble_quad(mesh,mat)
        n_cells=len(mesh.cells)
    else:
        mesh=make_mesh(**dict(g,element_type=kind))
        Kt,Ct=assemble_thermal(mesh,mat)
        Ku,_=assemble_elastic(mesh,mat,np.full(len(mesh.nodes),tref))
        M,G=assemble_inertia_coupling(mesh,mat)
        n_cells=len(mesh.elements)
    H=tref*G.T if not disable_feedback else 0*G.T
    n=len(mesh.nodes);nd=2*n;nodes=mesh.nodes
    bound=np.unique(np.concatenate(list(mesh.node_sets.values())))
    # Time-dependent exact BCs on all physical boundaries; prescribed as
    # mechanical x/y and theta DOFs. No face load: the Dirichlet reactions
    # account for the nonzero physical traction/heat flux.
    fixed=np.concatenate((2*bound,2*bound+1,nd+bound))
    fixed=np.sort(fixed)
    free=np.setdiff1d(np.arange(nd+n),fixed)
    beta=.25;gamma=.5;mcoef=1/(beta*dt*dt);vcoef=gamma/(beta*dt)
    A=bmat([[Ku+mcoef*M,-G],[vcoef*H,Kt+Ct/dt]],format='csc')
    factor=splu(A[free,:][:,free])
    af=A[free,:][:,fixed]
    u,initialT=exact_fields(nodes,0,t_ref=tref,curvature=curvature)
    theta=initialT-tref
    v=exact_initial_velocity(nodes)
    acc=2*0.02*curvature*np.column_stack((nodes[:,0]**2,nodes[:,0]*nodes[:,1])).ravel()
    for i in range(1,steps+1):
        t=i*dt
        f,q=_loads(mesh,kind,mat,t,curvature)
        up=u+dt*v+dt**2*(.5-beta)*acc
        vp=v+dt*(1-gamma)*acc
        rhs=np.r_[f+mcoef*M@up,q+Ct@theta/dt-H@(vp-vcoef*up)]
        ue,Te=exact_fields(nodes,t,t_ref=tref,curvature=curvature)
        known=np.r_[ue,Te-tref][fixed]
        out=np.zeros(nd+n);out[fixed]=known
        out[free]=factor.solve(rhs[free]-af@known)
        u_new=out[:nd];theta_new=out[nd:]
        acc_new=mcoef*(u_new-up)
        v=vp+gamma*dt*acc_new
        u,theta,acc=u_new,theta_new,acc_new
    u_exact,T_exact=exact_fields(nodes,end_time,t_ref=tref,curvature=curvature)
    dT=theta+tref-T_exact
    du=(u-u_exact).reshape(-1,2)
    return ManufacturedResult(kind,subdivisions,n,n_cells,steps,
                              float(np.sqrt(np.mean(dT**2))),
                              float(np.sqrt(np.mean(np.sum(du**2,axis=1)))),
                              float(np.max(abs(dT))),
                              float(np.max(np.linalg.norm(du,axis=1))),
                              float(np.ptp(T_exact)),
                              float(np.max(np.linalg.norm(u_exact.reshape(-1,2),axis=1))),
                              theta+tref,u,T_exact,u_exact)
