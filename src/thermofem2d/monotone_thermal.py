"""Experimental thermal-only discrete-maximum-principle reference.

For T3 and independent Q4, the conductivity matrix retains the original
spatial discretization unless a positive off-diagonal coupling would violate
monotonicity; such couplings are cancelled by a symmetric graph-Laplacian
correction. The heat-capacity matrix is row-sum lumped.

For T6, ordinary six-node quadratic shape functions have ZERO integral at
corner nodes. Consequently, direct row-sum lumping is singular and cannot
be used as a positivity-preserving scheme. A conforming 4x T3 subtriangle
space on the actual T6 nodes is used for BOTH capacity and conductivity.
It has a *different thermal trial space* and must not be represented as an
unchanged quadratic T6 temperature solution.

The graph correction is artificial diffusion. This is a first-order,
energy-conserving (row sums zero), positivity-preserving reference, NOT a
high-order flux-corrected transport scheme. It is restricted to isotropic
homogeneous heat conduction with no volumetric source or coupling feedback
that would invalidate a thermal-only maximum principle.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix, diags
from .element import plotting_triangles
from .fem import assemble_thermal
from .q4_reference import assemble_quad


@dataclass(frozen=True)
class MonotoneThermal:
    conductivity: object
    capacity: object
    initial_positive_offdiag_count: int
    added_graph_diffusion_L1: float
    capacity_total: float
    thermal_method: str


def _linear_subtriangle_heat_matrices(mesh, mat):
    """P1 assembly on four subtriangles per affine T6; no T6 mass lumping."""
    tri=plotting_triangles(mesh)
    xyz=mesh.nodes
    n=len(xyz)
    ii=[];jj=[];kv=[]
    cap=np.zeros(n,dtype=float)
    cond=float(mat['conductivity'])*float(mat['thickness'])
    volumetric_capacity=(float(mat['density'])*float(mat['heat_capacity'])*
                          float(mat['thickness']))
    if min(cond,volumetric_capacity)<=0:
        raise ValueError('positive thermal coefficients required')
    for nodes in tri:
        p=xyz[nodes]
        edge=np.stack((p[1]-p[0],p[2]-p[0]),axis=1)
        det=float(np.linalg.det(edge))
        if det<=1e-17:
            raise ValueError('degenerate T6 subtriangle')
        area=det/2
        dN=np.array([[-1.,-1.],[1.,0.],[0.,1.]]) @ np.linalg.inv(edge)
        stiff=cond*area*(dN@dN.T)
        cap[nodes]+=volumetric_capacity*area/3
        for a,i in enumerate(nodes):
            for b,j in enumerate(nodes):
                ii.append(int(i));jj.append(int(j));kv.append(float(stiff[a,b]))
    K=coo_matrix((kv,(ii,jj)),shape=(n,n)).tocsr()
    return K, cap


def _graph_monotonize(K):
    """Symmetric graph-Laplacian correction cancels positive off-diagonals."""
    K=K.tocsr()
    if K.shape[0]!=K.shape[1]:
        raise ValueError('square conductivity matrix expected')
    k=K.tocoo()
    selection=(k.row!=k.col)&(k.data>0)
    pi,pj,pv=k.row[selection],k.col[selection],k.data[selection]
    # Preserve strict symmetry to floating-point rounding.
    pos=coo_matrix((pv,(pi,pj)),shape=K.shape).tocsr()
    pos=(pos+pos.T)*0.5
    row=np.asarray(pos.sum(axis=1)).ravel()
    out=(K+diags(row)-pos).tocsr()
    out=(out+out.T)*0.5
    return out,int(np.count_nonzero(selection)),float(pv.sum())


def assemble_monotone_thermal(mesh, mat, family):
    """Return conservative monotone thermal matrices and traceable metadata.

    For Q4, the independent quadrilateral kernel is retained; for T3,
    the production linear-triangle kernel is retained. For T6, the thermal
    operator is replaced by an explicitly labelled P1-subtriangle operator.
    """
    if family=='T6':
        K,lump=_linear_subtriangle_heat_matrices(mesh,mat)
        method='P1 subtriangles on T6 nodes; row-lumped capacity'
    elif family=='T3':
        K,C=assemble_thermal(mesh,mat)
        lump=np.asarray(C.sum(axis=1)).ravel()
        method='T3 row-lumped capacity'
    elif family=='Q4':
        K,C,*_=assemble_quad(mesh,mat)
        lump=np.asarray(C.sum(axis=1)).ravel()
        method='Q4 row-lumped capacity'
    else:
        raise ValueError('family must be T3, Q4, or T6')
    if np.any(lump<=0) or not np.all(np.isfinite(lump)):
        raise ValueError('thermal mass lumping must be strictly positive')
    mono,npos,added=_graph_monotonize(K)
    rows=np.asarray(mono.sum(axis=1)).ravel()
    scale=max(1.,float(abs(mono.diagonal()).max()))
    if abs(rows).max()>1e-9*scale:
        raise ValueError('thermal operator has nonzero row sums')
    c=diags(lump,format='csr')
    return MonotoneThermal(mono,c,npos,added,float(lump.sum()),method)


def assemble_subtriangle_thermal_coupling(mesh, mat):
    """Consistent P2-displacement / subtriangle-P1-temperature coupling.

    The temperature interpolation is the same conforming P1 subtriangulation
    used in ``assemble_monotone_thermal`` for affine T6 geometry. The mechanical
    strain field retains the twelve-DOF quadratic T6 approximation. Reciprocal
    energy coupling uses the transpose of this matrix multiplied by T_ref.
    """
    from .fem import plane_stress_matrix
    from .element import shape_and_gradients, strain_displacement
    from scipy.sparse import coo_matrix

    if mesh.elements.shape[1]!=6:
        raise ValueError('P2/P1 coupling requires a six-node triangle mesh')
    dth=plane_stress_matrix(float(mat['E']),float(mat['nu'])) @ np.array(
        [float(mat['alpha']),float(mat['alpha']),0.]
    )
    thickness=float(mat['thickness'])
    n=len(mesh.nodes)
    indices_i=[];indices_j=[];vals=[]
    subtris=((0,3,5),(3,1,4),(5,4,2),(3,4,5))
    quad=np.array([[2./3.,1./6.,1./6.],
                   [1./6.,2./3.,1./6.],
                   [1./6.,1./6.,2./3.]])
    for conn in mesh.elements:
        xy=mesh.nodes[conn]
        parent=xy[:3]
        jac=np.stack((parent[1]-parent[0],parent[2]-parent[0]),axis=1)
        invjac=np.linalg.inv(jac)
        local=np.zeros((12,6),dtype=float)
        for local_nodes in subtris:
            p=xy[np.array(local_nodes)]
            area=abs(np.linalg.det(np.stack((p[1]-p[0],p[2]-p[0]),axis=1)))/2
            if area<=1e-17:raise ValueError('degenerate thermal subtriangle')
            for bary in quad:
                x=bary@p
                r,s=invjac@(x-parent[0])
                _,grad,_=shape_and_gradients(parent,'T6',float(r),float(s))
                Bt=strain_displacement(grad).T@dth
                local[:,local_nodes]+=thickness*area/3*np.outer(Bt,bary)
        displacement_dofs=np.ravel(np.column_stack((2*conn,2*conn+1)))
        for a,ia in enumerate(displacement_dofs):
            for b,ib in enumerate(conn):
                indices_i.append(int(ia));indices_j.append(int(ib));vals.append(float(local[a,b]))
    return coo_matrix((vals,(indices_i,indices_j)),shape=(2*n,n)).tocsr()
