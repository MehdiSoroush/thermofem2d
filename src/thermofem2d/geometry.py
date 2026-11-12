"""Conforming T3 mesh for a rectangular plate with one rectangular void.

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
    elements: np.ndarray  # (M,3) anticlockwise node IDs
    edges: dict[str, np.ndarray]  # boundary edges, IDs with shape (n,2)
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
              per_side=12, radial=12) -> Mesh:
    """Create a geometry-fitted mesh; coordinates in metres.

    Hole is axis-aligned; all boundaries follow exact straight edges. Meshing
    parameters must be integers. A small positive ligament is required.
    """
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
    return Mesh(nodes=coords,elements=tris,edges=edge_arrays,node_sets=node_sets)
