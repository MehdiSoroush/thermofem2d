"""T3 finite-element matrices and boundary operators.

Theory: small-strain, isotropic plane stress, Fourier heat conduction.
No thermo-mechanical feedback is included in this first release.
"""
from collections.abc import Mapping
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from .geometry import Mesh, BOUNDARIES


def _validate_material(material: Mapping) -> None:
    positive=("E","conductivity","density","heat_capacity","thickness")
    if any(not np.isfinite(material[k]) or material[k]<=0 for k in positive):
        raise ValueError("E, conductivity, density, heat_capacity, thickness must be positive")
    if not (-1 < material["nu"] < .5):
        raise ValueError("Poisson's ratio must lie in (-1, 0.5)")
    if not np.isfinite(material["alpha"]) or not np.isfinite(material["T_ref"]):
        raise ValueError("alpha and T_ref must be finite")


def plane_stress_matrix(E: float, nu: float) -> np.ndarray:
    return E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])


def _element_geometry(p: np.ndarray):
    x,y=p[:,0],p[:,1]
    twice=(x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0])
    if twice <= 0:
        raise ValueError("elements must be counter-clockwise and nondegenerate")
    grads=np.array([[y[1]-y[2],x[2]-x[1]],
                    [y[2]-y[0],x[0]-x[2]],
                    [y[0]-y[1],x[1]-x[0]]])/twice
    return twice/2,grads


def _assemble(ndof, rows, cols, values):
    return coo_matrix((values,(rows,cols)),shape=(ndof,ndof)).tocsr()


def assemble_thermal(mesh: Mesh, mat: Mapping):
    """Return conductivity K and consistent heat-capacity matrix C."""
    _validate_material(mat)
    rows=[]; cols=[]; kv=[]; cv=[]
    h=mat["thickness"]; cond=mat["conductivity"]
    rhoc=mat["density"]*mat["heat_capacity"]
    for tri in mesh.elements:
        area,g=_element_geometry(mesh.nodes[tri])
        Ke=cond*h*area*(g@g.T)
        Ce=rhoc*h*area/12*np.array([[2,1,1],[1,2,1],[1,1,2]])
        for a,i in enumerate(tri):
            for b,j in enumerate(tri):
                rows.append(int(i));cols.append(int(j));kv.append(Ke[a,b]);cv.append(Ce[a,b])
    n=len(mesh.nodes)
    return _assemble(n,rows,cols,kv),_assemble(n,rows,cols,cv)


def _boundary_conditions(mesh, specs, ndof, component=None):
    """Map boundary labels to nodal fixed DOFs, checking corner consistency."""
    fixed={}
    for label,value in specs.items():
        if label not in BOUNDARIES:
            raise ValueError(f"unknown boundary {label!r}")
        if not np.isfinite(value):
            raise ValueError("boundary values must be finite")
        for nid in mesh.node_sets[label]:
            dof=int(nid) if component is None else 2*int(nid)+component
            if dof in fixed and not np.isclose(fixed[dof],value,atol=1e-10,rtol=1e-10):
                raise ValueError(f"conflicting prescribed values at DOF {dof}")
            fixed[dof]=float(value)
    return fixed


def prescribed_temperature(mesh, specs):
    return _boundary_conditions(mesh,specs,len(mesh.nodes))


def prescribed_displacements(mesh, specs):
    fixed={}
    for component, axis in (("ux",0),("uy",1)):
        part=_boundary_conditions(mesh,specs.get(component,{}),2*len(mesh.nodes),axis)
        fixed.update(part)
    extras=set(specs)-{"ux","uy"}
    if extras:
        raise ValueError(f"unknown displacement components: {sorted(extras)}")
    return fixed


def solve_dirichlet(A:csr_matrix, rhs:np.ndarray, fixed:dict[int,float]):
    """Sparse essential-BC elimination, preserving physical prescribed values."""
    n=A.shape[0]
    if A.shape!=(n,n) or rhs.shape!=(n,):
        raise ValueError("incompatible matrix or right-hand side size")
    if not fixed:
        raise ValueError("at least one essential boundary constraint is required")
    ids=np.array(sorted(fixed),dtype=int)
    if ids.min()<0 or ids.max()>=n:
        raise ValueError("constrained DOF is outside the system")
    u=np.zeros(n)
    u[ids]=[fixed[k] for k in ids]
    free=np.setdiff1d(np.arange(n),ids)
    if free.size:
        u[free]=spsolve(A[free,:][:,free],rhs[free]-A[free,:][:,ids]@u[ids])
        if not np.all(np.isfinite(u[free])):
            raise ValueError("linear system singular or failed to converge; check constraints")
    return u


def boundary_load(mesh, loads, multiplier=1, components=1):
    """Integrate constant-in-space edge flux/traction with exact T3 line loads."""
    result=np.zeros(len(mesh.nodes)*components)
    for label,value in loads.items():
        if label not in BOUNDARIES:
            raise ValueError(f"unknown boundary {label!r}")
        val=np.asarray(value,dtype=float).reshape(-1)
        if len(val)!=components or not np.all(np.isfinite(val)):
            raise ValueError("incorrect boundary load components")
        for a,b in mesh.edges[label]:
            L=np.linalg.norm(mesh.nodes[b]-mesh.nodes[a])
            for nid in (a,b):
                result[nid*components:(nid+1)*components]+=multiplier*L/2*val
    return result


def assemble_elastic(mesh:Mesh, mat:Mapping, temperature:np.ndarray, tractions=None):
    """Assemble plane-stress stiffness and consistent thermal-expansion load."""
    _validate_material(mat)
    if temperature.shape!=(len(mesh.nodes),):
        raise ValueError("temperature vector shape mismatch")
    D=plane_stress_matrix(mat["E"],mat["nu"])
    rows=[];cols=[];vals=[]
    f=boundary_load(mesh, tractions or {},mat["thickness"],components=2)
    for tri in mesh.elements:
        area,g=_element_geometry(mesh.nodes[tri]); B=np.zeros((3,6))
        for n in range(3):
            dx,dy=g[n]; B[0,2*n]=dx;B[1,2*n+1]=dy
            B[2,2*n]=dy;B[2,2*n+1]=dx
        Ke=mat["thickness"]*area*(B.T@D@B)
        eps_th=mat["alpha"]*(temperature[tri].mean()-mat["T_ref"])*np.array([1,1,0])
        fe=mat["thickness"]*area*(B.T@D@eps_th)
        dof=np.ravel([[2*int(n),2*int(n)+1] for n in tri])
        f[dof]+=fe
        for a,i in enumerate(dof):
            for b,j in enumerate(dof):
                rows.append(i);cols.append(j);vals.append(Ke[a,b])
    return _assemble(2*len(mesh.nodes),rows,cols,vals),f


def stresses(mesh:Mesh,mat:Mapping,temperature:np.ndarray,displacement:np.ndarray):
    """Raw element-centroid stresses, not nodally smoothed quantities."""
    D=plane_stress_matrix(mat["E"],mat["nu"])
    sigma=np.empty((len(mesh.elements),3))
    for i,tri in enumerate(mesh.elements):
        _,g=_element_geometry(mesh.nodes[tri]);B=np.zeros((3,6))
        for n,(dx,dy) in enumerate(g):
            B[0,2*n]=dx;B[1,2*n+1]=dy
            B[2,2*n]=dy;B[2,2*n+1]=dx
        ue=displacement[np.ravel([[2*int(n),2*int(n)+1] for n in tri])]
        eth=mat["alpha"]*(temperature[tri].mean()-mat["T_ref"])*np.array([1,1,0])
        sigma[i]=D@(B@ue-eth)
    sx,sy,txy=sigma.T
    von_mises=np.sqrt(np.maximum(sx*sx-sx*sy+sy*sy+3*txy*txy,0))
    return sigma,von_mises


def element_heat_flux(mesh: Mesh, mat: Mapping, temperature: np.ndarray) -> np.ndarray:
    """Return elementwise Fourier heat flux q = -k grad(T), in W/m².

    The T3 temperature gradient is constant within each element.  The output
    represents flux per unit area, not an edge-integrated heat flow.  Input
    temperature is in kelvin and mesh geometry in metres.
    """
    _validate_material(mat)
    temperature = np.asarray(temperature, dtype=float)
    if temperature.shape != (len(mesh.nodes),) or not np.all(np.isfinite(temperature)):
        raise ValueError("expected one finite temperature per mesh node")
    flux = np.empty((len(mesh.elements), 2), dtype=float)
    for e, tri in enumerate(mesh.elements):
        _, gradients = _element_geometry(mesh.nodes[tri])
        flux[e] = -mat["conductivity"] * (temperature[tri] @ gradients)
    return flux
