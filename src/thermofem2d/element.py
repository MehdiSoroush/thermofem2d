"""Reference triangular T3/T6 Lagrange elements and degree-four Gauss rule.

Reference triangle: (r,s) = (0,0), (1,0), (0,1).
T6 order: three vertices followed by midsides 01, 12, 20.
All supported triangles are straight-sided (affine). Curved-edge quadratic
geometries require a nonconstant Jacobian and are outside this release.
"""
import numpy as np

# Six-point Dunavant degree-4 quadrature. Weights sum to one, so integrals
# on a triangle equal physical area times the weighted quadrature sum.
_A=0.816847572980459
_B=0.091576213509771
_C=0.108103018168070
_D=0.445948490915965
TRIANGLE_QUADRATURE = tuple(
    (p[1],p[2],w) for a,b,w in
    ((_A,_B,0.109951743655322),(_C,_D,0.223381589678011))
    for p in ((a,b,b),(b,a,b),(b,b,a))
)


def shape_and_gradients(vertices: np.ndarray, element_type: str, r: float, s: float):
    """Return N, physical dN/d(x,y), and positive triangle area.

    `vertices` is 3x2 for the oriented corner triangle; midside coordinates
    in the T6 mesh are required to coincide with edge midpoints.
    """
    x=np.asarray(vertices,dtype=float)
    if x.shape != (3,2):
        raise ValueError("expected 3 vertex coordinates")
    jac=np.stack((x[1]-x[0],x[2]-x[0]),axis=1)
    det=float(np.linalg.det(jac))
    if not np.isfinite(det) or det<=0:
        raise ValueError("triangle vertices must be counter-clockwise and nondegenerate")
    L=np.array([1-r-s,r,s],dtype=float)
    dL=np.array([[-1.,-1.],[1.,0.],[0.,1.]])@np.linalg.inv(jac)
    if element_type=="T3":
        return L,dL,0.5*det
    if element_type!="T6":
        raise ValueError("only T3 and T6 are supported")
    N=np.r_[L*(2*L-1),4*L[0]*L[1],4*L[1]*L[2],4*L[2]*L[0]]
    dN=np.vstack(((4*L-1)[:,None]*dL,
                  4*(L[0]*dL[1]+L[1]*dL[0]),
                  4*(L[1]*dL[2]+L[2]*dL[1]),
                  4*(L[2]*dL[0]+L[0]*dL[2])))
    return N,dN,0.5*det


def element_type_for_mesh(mesh):
    n=int(mesh.elements.shape[1])
    if n not in (3,6):
        raise ValueError("mesh supports only T3 or T6 triangles")
    return "T3" if n==3 else "T6"


def strain_displacement(gradients):
    """2D engineering-strain matrix [exx,eyy,gxy] for N displacement nodes."""
    g=np.asarray(gradients,dtype=float)
    B=np.zeros((3,2*len(g)))
    B[0,0::2]=g[:,0]
    B[1,1::2]=g[:,1]
    B[2,0::2]=g[:,1]
    B[2,1::2]=g[:,0]
    return B


def plotting_triangles(mesh):
    """Triangulate T6 visualization using all six nodes, without extra DOFs."""
    if element_type_for_mesh(mesh)=="T3":
        return mesh.elements
    e=mesh.elements
    return np.stack((e[:,[0,3,5]],e[:,[3,1,4]],
                     e[:,[5,4,2]],e[:,[3,4,5]]),axis=1).reshape(-1,3)
