"""Geometry and convergence checks for the circular-hole mesh utility."""

import numpy as np
import pytest

from thermofem2d.circular_geometry import (
    geometric_jacobian, hole_boundary_length, integrated_geometry_area,
    make_circular_mesh,
)
from thermofem2d.geometry import make_mesh


RADIUS = .25
CENTRE = (.12, -.07)


def mesh(kind, n, **kwargs):
    return make_circular_mesh(hole_radius=RADIUS, hole_cx=CENTRE[0],
                              hole_cy=CENTRE[1], per_side=n, radial=n,
                              element_type=kind, **kwargs)


@pytest.mark.parametrize("kind", ("T3", "T6"))
def test_circular_nodes_and_boundary_labels(kind):
    result = mesh(kind, 4)
    radius = np.linalg.norm(result.coordinates[result.boundary_nodes["hole"]] - CENTRE, axis=1)
    assert np.max(np.abs(radius-RADIUS)) < 1e-14
    assert len(result.boundary_edges["hole"]) == 16
    assert result.triangles.shape[1] == (3 if kind == "T3" else 6)
    assert set(result.boundary_edges) == {"hole", "bottom", "top", "left", "right"}
    assert all(len(result.boundary_edges[label]) == 4 for label in ("bottom", "top", "left", "right"))
    assert len(np.unique(result.boundary_edges["hole"])) == (16 if kind=="T3" else 32)


def test_quadratic_curved_midnodes_shared_by_adjacent_elements():
    result = mesh("T6", 4)
    for a, mid, b in result.boundary_edges["hole"]:
        midpoint = 0.5*(result.coordinates[a] + result.coordinates[b])
        assert np.linalg.norm(result.coordinates[mid]-midpoint) > 1e-4
        assert np.isclose(np.linalg.norm(result.coordinates[mid]-CENTRE), RADIUS)
        adjacent = [conn for conn in result.triangles
                    if a in conn[:3] and b in conn[:3]]
        assert len(adjacent) == 1
        assert mid in adjacent[0][3:]


@pytest.mark.parametrize("kind", ("T3", "T6"))
def test_jacobians_positive_for_eccentric_holes(kind):
    result = mesh(kind, 5)
    for conn in result.triangles:
        for r,s in ((0.,0.), (1.,0.), (0.,1.), (1/3,1/3), (.5,0.),(.5,.5),(0.,.5)):
            assert geometric_jacobian(result.coordinates,conn,r,s) > 0


def test_area_perimeter_convergence_and_order():
    actual_area = 4 - np.pi*RADIUS**2
    actual_length = 2*np.pi*RADIUS
    errors = {}
    for kind in ("T3", "T6"):
        errors[kind] = []
        for n in (3, 6, 10):
            result = mesh(kind,n)
            errors[kind].append((abs(integrated_geometry_area(result)-actual_area),
                                 abs(hole_boundary_length(result)-actual_length)))
        assert errors[kind][-1][0] < errors[kind][0][0]
        assert errors[kind][-1][1] < errors[kind][0][1]
    assert errors["T6"][-1][0] < errors["T3"][-1][0]/1000
    assert errors["T6"][-1][1] < errors["T3"][-1][1]/1000


def test_not_solver_ready_mesh():
    result = mesh("T6", 3)
    assert not hasattr(result, "nodes")
    assert not hasattr(result, "elements")


def test_backward_compatible_rectangular_geometry():
    result = make_mesh(per_side=4, radial=4, element_type="T6")
    assert result.elements.shape[1] == 6
    assert result.edges["hole"].shape == (16,3)


@pytest.mark.parametrize("options", [
    {"hole_radius": 0.}, {"hole_radius": -1.},
    {"hole_radius": np.nan}, {"hole_radius": np.inf},
    {"hole_radius": 1.0}, {"hole_cx": 0.9, "hole_radius": .2},
    {"per_side": 1}, {"radial": 0}, {"per_side": True},
    {"radial": 3.5}, {"element_type": "Q4"},
    {"hole_radius": "not a number"}, {"width": 0}, {"width": np.nan},
])
def test_bad_circular_geometry_rejected(options):
    with pytest.raises(ValueError):
        make_circular_mesh(**options)
