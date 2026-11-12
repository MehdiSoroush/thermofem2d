"""Conforming T3/T6 mesh for a rectangular plate with one rectangular void.

The annular mapping joins corresponding points on inner and outer boundaries.
It is structured topologically but yields an unstructured triangular FEM mesh.
The void does not contain elements and its boundary receives an explicit label.
"""
from dataclasses import dataclass
import numpy as np

BOUNDARIES = ("left", "right", "bottom", "top", "hole")


@dataclass(frozen=True)
class Mesh:
    nodes: np.ndarray  # (N,2) metre
    elements: np.ndarray  # (M,3) for T3 or (M,6) for T6
    edges: dict[str, np.ndarray]  # boundary edges: (n,2) T3 or (n,3) T6 [end,mid,end]
    node_sets: dict[str, np.ndarray]


def _corners(cx: float, cy: float, ax: float, ay: float) -> np.ndarray:
    return np.array([[cx-ax,cy-ay], [cx+ax,cy-ay],
                     [cx+ax,cy+ay], [cx-ax,cy+ay]], dtype=float)


def _perimeter_point(corners: np.ndarray, k: int, n: int) -> np.ndarray:
    side, i = divmod(k, n)
    r = i/n
    return (1-r)*corners[side] + r*corners[(side+1)%4]


def make_mesh(*, width=2.0, height=2.0, hole_width=0.24,
              hole_height=0.20, hole_cx=0.0, hole_cy=0.0,
              per_side=12, radial=12, element_type="T3") -> Mesh:
    """Create a geometry-fitted mesh; coordinates in metres.

    Hole is axis-aligned; all boundaries follow exact straight edges. Meshing
    parameters must be integers. A small positive ligament is required.
    """
    if element_type not in ("T3", "T6"):
        raise ValueError("element_type must be T3 or T6")
    if isinstance(per_side, bool) or not isinstance(per_side, int) or per_side<2:
        raise ValueError("per_side must be an integer >= 2")
    if isinstance(radial, bool) or not isinstance(radial, int) or radial<2:
        raise ValueError("radial must be an integer >= 2")
    vals=(width,height,hole_width,hole_height,hole_cx,hole_cy)
    if not np.all(np.isfinite(vals)):
        raise ValueError("geometry parameters must be finite")
    if width<=0 or height<=0 or hole_width<=0 or hole_height<=0:
        raise ValueError("dimensions must be positive")
    clearance_x=width/2 - (abs(hole_cx)+hole_width/2)
    clearance_y=height/2 - (abs(hole_cy)+hole_height/2)
    if min(clearance_x,clearance_y) <= 1e-9*max(width,height):
        raise ValueError("hole must lie strictly inside the plate")
    outer=_corners(0,0,width/2,height/2)
    inner=_corners(hole_cx,hole_cy,hole_width/2,hole_height/2)
    count=4*per_side
    coords=np.empty(((radial+1)*count, 2),dtype=float)
    for j in range(radial+1):
        t=j/radial
        for k in range(count):
            coords[j*count+k]=(1-t)*_perimeter_point(inner,k,per_side)+t*_perimeter_point(outer,k,per_side)
    cells=[]
    for j in range(radial):
        for k in range(count):
            kn=(k+1)%count
            a=j*count+k; b=j*count+kn
            c=(j+1)*count+kn; d=(j+1)*count+k
            cells.extend(((a,b,c),(a,c,d)))
    tris=np.array(cells,dtype=int)
    p=coords[tris]
    signed=(p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1])-(p[:,1,1]-p[:,0,1])*(p[:,2,0]-p[:,0,0])
    if np.any(np.isclose(signed,0,atol=1e-14*width*height)):
        raise ValueError("mesh contains degenerate triangles")
    flip=signed<0
    tris[flip,1],tris[flip,2]=tris[flip,2].copy(),tris[flip,1].copy()
    edges={key:[] for key in BOUNDARIES}
    for k in range(count):
        side=k//per_side
        kn=(k+1)%count
        edges[("bottom","right","top","left")[side]].append((radial*count+k,radial*count+kn))
        edges["hole"].append((k,kn))
    edge_arrays={key:np.asarray(v,dtype=int).reshape(-1,2) for key,v in edges.items()}
    node_sets={key:np.unique(e) for key,e in edge_arrays.items()}
    basic=Mesh(nodes=coords,elements=tris,edges=edge_arrays,node_sets=node_sets)
    return basic if element_type == "T3" else upgrade_to_t6(basic)


def upgrade_to_t6(mesh: Mesh) -> Mesh:
    """Raise a straight-sided conforming T3 mesh to quadratic T6.

    Local order [v0,v1,v2,m01,m12,m20]. Global midside nodes are
    shared by neighboring elements, including boundary edges. This is an
    *affine* T6 mapping (midsides at edge midpoints), not a curved-edge CAD
    approximation. Conversion preserves the original polygonal geometry.
    """
    if mesh.elements.ndim != 2 or mesh.elements.shape[1] != 3:
        raise ValueError("upgrade_to_t6 expects a T3 mesh")
    pts=[p.copy() for p in mesh.nodes]
    mids={}
    def midpoint(a,b):
        key=(min(int(a),int(b)),max(int(a),int(b)))
        if key not in mids:
            mids[key]=len(pts)
            pts.append(0.5*(mesh.nodes[key[0]]+mesh.nodes[key[1]]))
        return mids[key]
    e6=[]
    for a,b,c in mesh.elements:
        e6.append((a,b,c,midpoint(a,b),midpoint(b,c),midpoint(c,a)))
    boundary={}
    for label,edges in mesh.edges.items():
        boundary[label]=np.array([(a,midpoint(a,b),b) for a,b in edges],dtype=int).reshape(-1,3)
    return Mesh(nodes=np.asarray(pts,dtype=float), elements=np.asarray(e6,dtype=int),
                edges=boundary, node_sets={label:np.unique(e) for label,e in boundary.items()})
