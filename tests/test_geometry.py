import numpy as np
import pytest
from thermofem2d.geometry import make_mesh
from thermofem2d.fem import _element_geometry


def test_mesh_closed_hole():
    m=make_mesh(per_side=7,radial=5)
    assert len(m.nodes)==6*4*7
    assert len(m.elements)==2*5*4*7
    assert len(m.edges["hole"])==28
    assert set(m.edges["hole"].reshape(-1))==set(m.node_sets["hole"])
    assert all(_element_geometry(m.nodes[tri])[0]>0 for tri in m.elements)


def test_offset_hole_remains_empty():
    m=make_mesh(hole_cx=.3,hole_cy=-.18,hole_width=.4,hole_height=.2)
    centers=m.nodes[m.elements].mean(axis=1)
    inside=(abs(centers[:,0]-.3)<.20)&(abs(centers[:,1]+.18)<.10)
    assert not inside.any()
    hole=m.nodes[m.node_sets["hole"]]
    assert min(hole[:,0])==pytest.approx(.1)
    assert max(hole[:,0])==pytest.approx(.5)


@pytest.mark.parametrize("kwargs",[
    dict(hole_cx=.95),dict(hole_width=2),dict(hole_height=0),
    dict(per_side=1),dict(radial=0),dict(width=-1),dict(hole_cx=float("nan"))])
def test_invalid_geometry_rejected(kwargs):
    with pytest.raises(ValueError):make_mesh(**kwargs)
