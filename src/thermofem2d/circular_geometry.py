"""Geometric circular-hole T3/T6 meshes and independent area diagnostics.

This module is a geometry-only utility. The production solver currently
supports straight-sided elements and rectangular openings. Its affine
matrices must not be applied to the nonaffine circular geometry.

Coordinates are expressed in metres. T6 hole-edge nodes lie on the specified
circle, but quadratic Lagrange edges approximate rather than represent exact
circular arcs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import Mesh, _corners, _perimeter_point, upgrade_to_t6


@dataclass(frozen=True)
class CircularGeometry:
    """Geometric data deliberately separate from solver-ready ``Mesh``.

    This type cannot be passed directly to the production FEM assembler.
    Boundary-edge order for T6 is [endpoint, midside, endpoint].
    """
    coordinates: np.ndarray
    triangles: np.ndarray
    boundary_edges: dict[str, np.ndarray]
    boundary_nodes: dict[str, np.ndarray]



def _validate_dimensions(width, height, radius, hole_cx, hole_cy, per_side, radial, element_type):
    if element_type not in ("T3", "T6"):
        raise ValueError("element_type must be T3 or T6")
    for label, value in (("per_side", per_side), ("radial", radial)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 2:
            raise ValueError(f"{label} must be an integer >= 2")
    dims = (width, height, radius, hole_cx, hole_cy)
    if not all(np.isscalar(x) and not isinstance(x, (str, bytes, bool)) for x in dims):
        raise ValueError("geometry dimensions must be finite real numbers")
    try:
        vals = np.asarray(dims, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("geometry dimensions must be finite real numbers") from exc
    if not np.all(np.isfinite(vals)):
        raise ValueError("geometry dimensions must be finite")
    width, height, radius, hole_cx, hole_cy = vals
    if min(width, height, radius) <= 0:
        raise ValueError("plate dimensions and radius must be positive")
    clearance_x = width / 2 - abs(hole_cx) - radius
    clearance_y = height / 2 - abs(hole_cy) - radius
    if min(clearance_x, clearance_y) <= 1e-9 * max(width, height):
        raise ValueError("circular hole must lie strictly inside the plate")
    return width, height, radius, hole_cx, hole_cy


def _reference_shape_gradients(kind, r, s):
    """Natural-coordinate derivatives of the T3/T6 Lagrange shape functions."""
    L = np.array([1.0-r-s, r, s], dtype=float)
    dL = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    if kind == "T3":
        return dL
    return np.vstack(((4*L-1)[:, None] * dL,
                      4*(L[0]*dL[1] + L[1]*dL[0]),
                      4*(L[1]*dL[2] + L[2]*dL[1]),
                      4*(L[2]*dL[0] + L[0]*dL[2])))


def geometric_jacobian(nodes, connectivity, r, s):
    """Compute det(dx/dr,dx/ds) from all element geometry nodes.

    This is a geometry diagnostic, not a substitute for physical FEM assembly.
    """
    points = np.asarray(nodes, dtype=float)[np.asarray(connectivity, dtype=int)]
    if points.shape[0] not in (3, 6) or points.shape[1] != 2:
        raise ValueError("triangular element connectivity must have 3 or 6 nodes")
    if not np.all(np.isfinite(points)):
        raise ValueError("element coordinates must be finite")
    kind = "T3" if len(points) == 3 else "T6"
    return float(np.linalg.det(points.T @ _reference_shape_gradients(kind, r, s)))


def _sample_points():
    # A 5x5 mapped quadrature set plus corners and edge midpoints samples
    # the nonaffine Jacobian at interior and boundary locations.
    gauss, weights = np.polynomial.legendre.leggauss(5)
    x = 0.5 * (gauss + 1)
    w = 0.5 * weights
    values = [(float(u), float((1-u)*v), float(2*wu*wv*(1-u)))
              for u, wu in zip(x, w) for v, wv in zip(x, w)]
    samples = values + [(r, s, 0.0) for r, s in
                        ((0., 0.), (1., 0.), (0., 1.), (.5, 0.), (.5, .5), (0., .5))]
    return samples


_QUADRATURE_AND_SAMPLES = _sample_points()


def integrated_geometry_area(mesh):
    """Numerically integrate the geometric map over all triangles (m²).

    Weights are normalized so that physical area = sum(det(J) * weight / 2).
    """
    total = 0.0
    for conn in mesh.triangles:
        for r, s, w in _QUADRATURE_AND_SAMPLES:
            if w:
                total += 0.5 * w * geometric_jacobian(mesh.coordinates, conn, r, s)
    return float(total)


def hole_boundary_length(mesh):
    """Integrate the piecewise linear/quadratic hole boundary arc (m)."""
    xi, weights = np.polynomial.legendre.leggauss(5)
    length = 0.0
    for conn in mesh.boundary_edges["hole"]:
        p = mesh.coordinates[conn]
        if len(conn) == 2:
            length += float(np.linalg.norm(p[1] - p[0]))
        else:
            for z, w in zip(xi, weights):
                dx = (z-0.5)*p[0] - 2*z*p[1] + (z+0.5)*p[2]
                length += w * float(np.linalg.norm(dx))
    return float(length)


def make_circular_mesh(*, width=2.0, height=2.0, hole_radius=0.25,
                       hole_cx=0.0, hole_cy=0.0, per_side=8, radial=8,
                       element_type="T3") -> CircularGeometry:
    """Create a square/rectangular plate mesh with one circular through-hole.

    The returned geometry has a distinct type rather than the production solver
    ``Mesh`` type: no thermal, mechanical, or coupled solution is provided.
    """
    width, height, radius, cx, cy = _validate_dimensions(
        width, height, hole_radius, hole_cx, hole_cy, per_side, radial, element_type)
    count = 4 * int(per_side)
    centre = np.array([cx, cy])
    outer = _corners(0.0, 0.0, width/2, height/2)
    points = np.empty(((int(radial)+1)*count, 2), dtype=float)
    for j in range(int(radial)+1):
        t = j/int(radial)
        for k in range(count):
            phi = -3*np.pi/4 + 2*np.pi*k/count
            inner = centre + radius * np.array([np.cos(phi), np.sin(phi)])
            outside = _perimeter_point(outer, k, int(per_side))
            points[j*count+k] = (1-t)*inner+t*outside
    cells = []
    for j in range(int(radial)):
        for k in range(count):
            kn = (k+1) % count
            a, b = j*count+k, j*count+kn
            c, d = (j+1)*count+kn, (j+1)*count+k
            cells.extend(((a,b,c), (a,c,d)))
    tris = np.asarray(cells, dtype=int)
    p = points[tris]
    signed = ((p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1])-
              (p[:,1,1]-p[:,0,1])*(p[:,2,0]-p[:,0,0]))
    if np.any(signed == 0) or not np.all(np.isfinite(signed)):
        raise ValueError("degenerate circular-hole mesh")
    flip = signed < 0
    tris[flip, 1], tris[flip, 2] = tris[flip, 2].copy(), tris[flip, 1].copy()
    boundary = {key: [] for key in ("left", "right", "top", "bottom", "hole")}
    for k in range(count):
        kn = (k+1) % count
        label = ("bottom", "right", "top", "left")[k//int(per_side)]
        boundary[label].append((int(radial)*count+k, int(radial)*count+kn))
        boundary["hole"].append((k,kn))
    edge_arrays = {k: np.asarray(v, dtype=int) for k,v in boundary.items()}
    mesh = Mesh(points, tris, edge_arrays,
                {k: np.unique(v) for k,v in edge_arrays.items()})
    if element_type == "T6":
        mesh = upgrade_to_t6(mesh)
        nodes = mesh.nodes.copy()
        for a, mid, b in mesh.edges["hole"]:
            direction = nodes[a] + nodes[b] - 2*centre
            norm = np.linalg.norm(direction)
            if norm <= 1e-12 * radius:
                raise ValueError("circular edge subtends an invalid angle")
            nodes[mid] = centre + radius * direction/norm
        mesh = Mesh(nodes, mesh.elements, mesh.edges, mesh.node_sets)
    # A dedicated type prevents accidental solver dispatch; sampled Jacobian
    # positivity is necessary, but not a formal global validity proof.
    geometry = CircularGeometry(mesh.nodes, mesh.elements, mesh.edges, mesh.node_sets)
    scale = max(width, height)
    for conn in geometry.triangles:
        for r, s, _ in _QUADRATURE_AND_SAMPLES:
            det = geometric_jacobian(geometry.coordinates, conn, r, s)
            if not np.isfinite(det) or det <= 1e-15 * scale**2:
                raise ValueError("invalid isoparametric triangle geometry")
    return geometry
