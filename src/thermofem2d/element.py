"""Reference triangular Lagrange elements and isoparametric integration.

Reference triangle: (r,s) = (0,0), (1,0), (0,1).
T6 order: three vertices followed by midsides 01, 12, 20.
T6 elements admit affine or non-affine (curved-edge) coordinates, using
the full six-node isoparametric Jacobian at every integration point.
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


# A 5x5 Gauss-Legendre tensor-product rule maps the unit square to the
# reference triangle: r=u, s=(1-u)v, dA=(1-u)du dv. Weights are normalized
# to sum to one, matching TRIANGLE_QUADRATURE conventions. The higher-order
# rule resolves curved Jacobians and geometric mass integration.
def _curved_triangle_quadrature():
    abscissae, weights=np.polynomial.legendre.leggauss(5)
    u=0.5*(abscissae+1)
    w=0.5*weights
    return tuple((float(a),float((1-a)*b),float(2*wa*wb*(1-a)))
                 for a,wa in zip(u,w) for b,wb in zip(u,w))


CURVED_TRIANGLE_QUADRATURE=_curved_triangle_quadrature()


def has_curved_elements(mesh, *, atol=1e-12):
    """Detect quadratic geometry that differs from straight-edge midpoints."""
    if mesh.elements.shape[1] != 6:
        return False
    e=mesh.elements
    corners=mesh.nodes[e[:,[0,1,1,2,2,0]]].reshape(-1,3,2,2)
    expected=corners.mean(axis=2)
    actual=mesh.nodes[e[:,3:]]
    span=max(1.0,float(np.max(np.ptp(mesh.nodes,axis=0))))
    return bool(np.max(np.abs(expected-actual))>atol*span)


def integration_rule(mesh):
    """Retain the legacy affine rule; use higher-order geometry quadrature."""
    return CURVED_TRIANGLE_QUADRATURE if has_curved_elements(mesh) else TRIANGLE_QUADRATURE


def shape_and_gradients(vertices: np.ndarray, element_type: str, r: float, s: float):
    """Return shape values, physical gradients, and local Jacobian area factor.

    For T6, accepts either three vertices (legacy affine geometry) or all six
    geometry nodes (isoparametric). The third value equals det(J)/2 and can
    vary across a curved element. Quadrature weights are normalized to one.
    """
    x=np.asarray(vertices,dtype=float)
    if element_type not in ("T3","T6"):
        raise ValueError("only T3 and T6 are supported")
    if x.shape!=(3,2) and not (element_type=="T6" and x.shape==(6,2)):
        raise ValueError("T3 requires 3 coordinates; T6 accepts 3 or 6")
    L=np.array([1-r-s,r,s],dtype=float)
    dL=np.array([[-1.,-1.],[1.,0.],[0.,1.]])
    if element_type=="T3":
        N=L
        dN_ref=dL
    else:
        N=np.r_[L*(2*L-1),4*L[0]*L[1],4*L[1]*L[2],4*L[2]*L[0]]
        dN_ref=np.vstack(((4*L-1)[:,None]*dL,
                          4*(L[0]*dL[1]+L[1]*dL[0]),
                          4*(L[1]*dL[2]+L[2]*dL[1]),
                          4*(L[2]*dL[0]+L[0]*dL[2])))
    # The T6 field has six basis functions even when supplied with three
    # corner coordinates. Only the geometric Jacobian then remains affine.
    geom_grad=dN_ref if len(x)==len(N) else dL
    jac=x.T@geom_grad
    det=float(np.linalg.det(jac))
    scale=max(1.0,float(np.max(np.abs(x))))
    if not np.isfinite(det) or det <= 1e-14*scale*scale:
        raise ValueError("nonpositive or degenerate isoparametric Jacobian")
    dN=dN_ref@np.linalg.inv(jac)
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
