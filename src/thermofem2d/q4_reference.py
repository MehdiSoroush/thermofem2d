"""Independent Q4 isoparametric FEM reference for the same *model class*.

Independently constructs a conforming ring-quad mesh, 2x2 Gauss shape functions,
thermal capacity/conductivity, plane-stress elasticity, consistent inertia,
reciprocal coupling, and boundary loads. Does not import production geometry,
element, FEM, or coupled assembly/stepper code. This is a second spatial
implementation, not an external solver or proof of continuum-level accuracy.

Supports the v0.5 subset: one rectangular through-hole, homogeneous isotropic
constant properties, constant scalar thermal and vector mechanical boundary
loads, constant Dirichlet data, zero initial displacement/velocity, and scalar
or affine temperature initial conditions. Positive imposed heat flux points
*into* the domain (the config convention of the original solver).
"""
from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np
from scipy.sparse import coo_matrix, bmat
from scipy.sparse.linalg import spsolve, splu


BOUNDARY_NAMES = ("bottom", "right", "top", "left", "hole")
_GAUSS = (-1 / np.sqrt(3), 1 / np.sqrt(3))
_SIGNS = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])


@dataclass
class QuadMesh:
    nodes: np.ndarray
    cells: np.ndarray
    edges: dict[str, np.ndarray]
    node_sets: dict[str, np.ndarray]


def quad_mesh(geometry: Mapping) -> QuadMesh:
    """Create Q4 cells independently of the production triangular mesh builder."""
    g=geometry
    width, height = float(g["width"]), float(g["height"])
    hw, hh = float(g["hole_width"]), float(g["hole_height"])
    cx, cy = float(g["hole_cx"]), float(g["hole_cy"])
    ns, nr = g["per_side"], g["radial"]
    if any(isinstance(z,bool) or not isinstance(z,int) or z<2 for z in (ns,nr)):
        raise ValueError("per_side and radial must be integers >=2")
    if not np.all(np.isfinite([width,height,hw,hh,cx,cy])) or min(width,height,hw,hh)<=0:
        raise ValueError("invalid geometry")
    if abs(cx)+hw/2>=width/2 or abs(cy)+hh/2>=height/2:
        raise ValueError("hole outside plate")
    outer=np.array([[-width/2,-height/2],[width/2,-height/2],
                    [width/2,height/2],[-width/2,height/2]])
    inner=np.array([[cx-hw/2,cy-hh/2],[cx+hw/2,cy-hh/2],
                    [cx+hw/2,cy+hh/2],[cx-hw/2,cy+hh/2]])
    count=4*ns
    nodes=[]
    for j in range(nr+1):
        t=j/nr
        for side in range(4):
            for k in range(ns):
                lam=k/ns
                inside=inner[side]*(1-lam)+inner[(side+1)%4]*lam
                outside=outer[side]*(1-lam)+outer[(side+1)%4]*lam
                nodes.append((1-t)*inside+t*outside)
    nodes=np.asarray(nodes)
    cells=[]
    for j in range(nr):
        for k in range(count):
            kn=(k+1)%count
            cells.append((j*count+k,(j+1)*count+k,(j+1)*count+kn,j*count+kn))
    edges={label:[] for label in BOUNDARY_NAMES}
    for k in range(count):
        kn=(k+1)%count
        edges[BOUNDARY_NAMES[k//ns]].append((nr*count+k,nr*count+kn))
        edges["hole"].append((k,kn))
    edges={label:np.asarray(conn,dtype=int) for label,conn in edges.items()}
    return QuadMesh(nodes,np.asarray(cells,dtype=int),edges,
                    {label:np.unique(conn) for label,conn in edges.items()})


def q4_shape(xi: float, eta: float):
    """N and gradients w.r.t. natural coordinates, independent of triangles."""
    x,y=_SIGNS.T
    N=.25*(1+x*xi)*(1+y*eta)
    dn=.25*np.column_stack((x*(1+y*eta),y*(1+x*xi)))
    return N,dn


def _matrix(rows,cols,values,n):
    return coo_matrix((values,(rows,cols)),shape=(n,n)).tocsr()


def assemble_quad(mesh:QuadMesh,mat:Mapping):
    """Return Ktheta,Ctheta,Ku,M,G from independent isoparametric Q4 assembly."""
    E,nu=float(mat["E"]),float(mat["nu"])
    if not 0<E or not -1<nu<.5:
        raise ValueError("unsupported elastic constants")
    k,rho,c,h,alpha=(float(mat[t]) for t in
                    ("conductivity","density","heat_capacity","thickness","alpha"))
    if min(k,rho,c,h)<=0:
        raise ValueError("thermal/mass parameters must be positive")
    D=E/(1-nu*nu)*np.array([[1.,nu,0],[nu,1.,0],[0,0,(1-nu)/2]])
    dth=D @ np.array([alpha,alpha,0.])
    n=len(mesh.nodes)
    tr,tc,tv,cv=[],[],[],[]
    mr,mc,mv,kv=[],[],[],[]
    gr,gc,gv=[],[],[]
    for conn in mesh.cells:
        xy=mesh.nodes[conn]
        Kth=np.zeros((4,4)); Cth=np.zeros((4,4))
        Ku=np.zeros((8,8)); M=np.zeros((8,8));G=np.zeros((8,4))
        for xi in _GAUSS:
            for eta in _GAUSS:
                N,ref_grad=q4_shape(xi,eta)
                jac=xy.T@ref_grad
                det=float(np.linalg.det(jac))
                if not np.isfinite(det) or det<=0:
                    raise ValueError("negative/inverted Q4 Jacobian")
                grad=ref_grad @ np.linalg.inv(jac)
                B=np.zeros((3,8))
                B[0,0::2]=grad[:,0]
                B[1,1::2]=grad[:,1]
                B[2,0::2]=grad[:,1]
                B[2,1::2]=grad[:,0]
                fac=det*h
                Kth+=fac*k*(grad@grad.T)
                Cth+=fac*rho*c*np.outer(N,N)
                Ku+=fac*B.T@D@B
                scalar_mass=fac*rho*np.outer(N,N)
                M[0::2,0::2]+=scalar_mass
                M[1::2,1::2]+=scalar_mass
                G+=fac*np.outer(B.T@dth,N)
        dofs=np.ravel(np.column_stack((2*conn,2*conn+1)))
        for a,i in enumerate(conn):
            for b,j in enumerate(conn):
                tr.append(int(i));tc.append(int(j))
                tv.append(Kth[a,b]);cv.append(Cth[a,b])
        for a,i in enumerate(dofs):
            for b,j in enumerate(dofs):
                mr.append(int(i));mc.append(int(j))
                kv.append(Ku[a,b]);mv.append(M[a,b])
            for b,j in enumerate(conn):
                gr.append(int(i));gc.append(int(j));gv.append(G[a,b])
    Kth=_matrix(tr,tc,tv,n);Cth=_matrix(tr,tc,cv,n)
    Ku=_matrix(mr,mc,kv,2*n);M=_matrix(mr,mc,mv,2*n)
    G=coo_matrix((gv,(gr,gc)),shape=(2*n,n)).tocsr()
    return Kth,Cth,Ku,M,G


def _fixed(mesh, specs, component=None):
    vals={}
    for label,value in specs.items():
        if label not in mesh.node_sets:
            raise ValueError("unknown boundary label")
        value=float(value)
        if not np.isfinite(value):
            raise ValueError("nonfinite prescribed value")
        for node in mesh.node_sets[label]:
            idx=int(node) if component is None else 2*int(node)+component
            if idx in vals and not np.isclose(vals[idx],value,rtol=1e-11,atol=1e-11):
                raise ValueError("inconsistent corner constraints")
            vals[idx]=value
    return vals


def _edge_load(mesh, specs, scale, components):
    vec=np.zeros(len(mesh.nodes)*components)
    for label,value in specs.items():
        if label not in mesh.edges:
            raise ValueError("unknown boundary label")
        f=np.asarray(value,dtype=float).reshape(-1)
        if len(f)!=components or not np.all(np.isfinite(f)):
            raise ValueError("incorrect boundary load")
        for e in mesh.edges[label]:
            p=mesh.nodes[e]
            length=float(np.linalg.norm(p[1]-p[0]))
            for node in e:
                vec[components*node:components*(node+1)]+=0.5*length*scale*f
    return vec


def _initial_temperature(mesh,analysis,mat):
    v=analysis.get("initial_temperature",mat["T_ref"])
    if isinstance(v,dict):
        if v.get("type")!="affine":
            raise ValueError("reference initial field requires affine specification")
        return np.array(v["base"]+v["gx"]*mesh.nodes[:,0]+v["gy"]*mesh.nodes[:,1],dtype=float)
    return np.full(len(mesh.nodes),float(v))


@dataclass
class QuadReference:
    mesh:QuadMesh
    temperature:np.ndarray
    displacement:np.ndarray
    times:np.ndarray
    temperature_history:np.ndarray
    displacement_history:np.ndarray
    feedback:bool


def _reduced_dirichlet(A,b,prescribed):
    out=np.zeros(A.shape[0]); ids=np.asarray(sorted(prescribed),dtype=int)
    free=np.setdiff1d(np.arange(len(out)),ids)
    if len(ids): out[ids]=[prescribed[int(x)] for x in ids]
    if len(free):
        out[free]=spsolve(A[free,:][:,free].tocsc(),b[free]-A[free,:][:,ids]@out[ids])
    if not np.all(np.isfinite(out)):
        raise ValueError("reference system is nonfinite")
    return out


def solve_q4_coupled(config:Mapping)->QuadReference:
    """Solve a second FEM spatial discretization of coupled transient physics."""
    a,g,mat=config["analysis"],config["geometry"],config["material"]
    if a.get("thermoelastic_mode")!="coupled_dynamic" or a.get("thermal_mode")!="transient":
        raise ValueError("Q4 coupled reference requires coupled_dynamic transient")
    if a.get("initial_displacement",0)!=0 or a.get("initial_velocity",0)!=0:
        raise ValueError("only zero initial displacement/velocity supported")
    t_ref=float(mat["T_ref"])
    if t_ref<=0: raise ValueError("T_ref must be absolute temperature")
    dt=float(a["time_step"]);steps=a["steps"]
    if not np.isfinite(dt) or dt<=0 or isinstance(steps,bool) or not isinstance(steps,int) or not 1<=steps<=10000:
        raise ValueError("invalid time grid")
    mesh=quad_mesh(g);N=len(mesh.nodes);nu=2*N
    Kt,Ct,Ku,M,G=assemble_quad(mesh,mat)
    feedback=a.get("thermoelastic_feedback",True)
    if not isinstance(feedback,bool):raise ValueError("feedback must be bool")
    H=t_ref*G.T if feedback else 0*G.T
    f=_edge_load(mesh,a.get("traction_bc",{}),mat["thickness"],2)
    q=_edge_load(mesh,a.get("inward_heat_flux",{}),mat["thickness"],1)
    Q=float(a.get("volumetric_heat_source",0.))
    if Q:
        for conn in mesh.cells:
            p=mesh.nodes[conn]
            for xi in _GAUSS:
                for eta in _GAUSS:
                    Nq,dn=q4_shape(xi,eta)
                    q[conn]+=Q*mat["thickness"]*np.linalg.det(p.T@dn)*Nq
    fu={}
    for name,component in (("ux",0),("uy",1)):
        fu.update(_fixed(mesh,a.get("displacement_bc",{}).get(name,{}),component))
    if not fu: raise ValueError("mechanical supports required")
    ft=_fixed(mesh,a.get("temperature_bc",{}))
    fixed=dict(fu);fixed.update({nu+i:v-t_ref for i,v in ft.items()})
    u=np.zeros(nu);v=np.zeros(nu)
    for idx,val in fu.items():u[idx]=val
    th=_initial_temperature(mesh,a,mat)-t_ref
    for idx,val in ft.items(): th[idx]=val-t_ref
    acc=np.zeros(nu)
    freeU=np.setdiff1d(np.arange(nu),list(fu))
    if len(freeU):
        acc[freeU]=spsolve(M[freeU,:][:,freeU].tocsc(),(f+G@th-Ku@u)[freeU])
    beta=.25;gamma=.5
    mcoef=1/(beta*dt**2);vcoef=gamma/(beta*dt)
    L=bmat([[Ku+mcoef*M,-G],[vcoef*H,Kt+Ct/dt]],format="csc")
    ids=np.asarray(sorted(fixed),dtype=int)
    free=np.setdiff1d(np.arange(nu+N),ids)
    vals=np.asarray([fixed[int(k)] for k in ids])
    factor=splu(L[free,:][:,free]) if len(free) else None
    fc=L[free,:][:,ids]
    ts=np.arange(steps+1)*dt
    t_hist=np.zeros((steps+1,N));u_hist=np.zeros((steps+1,nu))
    t_hist[0]=th+t_ref;u_hist[0]=u
    for k in range(1,steps+1):
        up=u+dt*v+dt**2*(.5-beta)*acc
        vp=v+dt*(1-gamma)*acc
        b=np.r_[f+mcoef*M@up,q+Ct@th/dt-H@(vp-vcoef*up)]
        x=np.zeros(nu+N);x[ids]=vals
        if len(free):x[free]=factor.solve(b[free]-fc@vals)
        nextu=x[:nu]; th=x[nu:]
        acc=mcoef*(nextu-up)
        v=vp+gamma*dt*acc;u=nextu
        t_hist[k]=th+t_ref;u_hist[k]=u
    return QuadReference(mesh,th+t_ref,u,ts,t_hist,u_hist,feedback)


def q4_patch_tests(mesh:QuadMesh,mat:Mapping):
    """Analytic affine temperature and displacement patch tests on quadrilaterals.

    For conductivity without sources, prescribed affine temperatures at *all*
    physical boundary nodes must reproduce a harmonic linear field. For
    constant-strain elasticity, affine displacements on all boundary nodes
    induce zero unbalanced force at interior nodes.
    """
    Kt,_,Ku,_,_=assemble_quad(mesh,mat)
    x,y=mesh.nodes.T;n=len(x)
    expectedT=310+8*x-3*y
    support=np.unique(np.concatenate(list(mesh.node_sets.values())))
    solvedT=_reduced_dirichlet(Kt,np.zeros(n),{int(i):expectedT[i] for i in support})
    ux=2e-5*x-1e-5*y; uy=-0.8e-5*x+1.3e-5*y
    expectedU=np.column_stack((ux,uy)).ravel()
    fixed={2*int(i)+j:expectedU[2*i+j] for i in support for j in range(2)}
    solvedU=_reduced_dirichlet(Ku,np.zeros(2*n),fixed)
    return float(np.max(np.abs(solvedT-expectedT))),float(np.max(np.abs(solvedU-expectedU)))
