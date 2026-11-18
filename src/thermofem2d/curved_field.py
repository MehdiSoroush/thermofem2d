"""Explicit adapter from verified circular geometry to isoparametric FEM assembly.

The circular-geometry builder remains independent of the production model
configuration. This adapter permits direct assembly of thermal and mechanical
matrices on a validated circular-hole T3/T6 mesh, without exposing curved
geometry to the normal CLI or coupled time-stepper model dispatch.
"""
from __future__ import annotations

import numpy as np

from .circular_geometry import CircularGeometry, geometric_jacobian
from .element import CURVED_TRIANGLE_QUADRATURE
from .geometry import Mesh


def circular_field_mesh(geometry: CircularGeometry) -> Mesh:
    """Validate and copy circular-hole geometry into the FEM assembly format.

    The T6 edge convention is [end0, midside, end1] and the T6 element
    convention is [vertex0, vertex1, vertex2, mid01, mid12, mid20]. The
    numerical Jacobian test checks every curved element at integration points
    and all three reference corners; it is not a global validity certificate.
    """
    if not isinstance(geometry, CircularGeometry):
        raise TypeError("expected a CircularGeometry instance")
    nodes = np.array(geometry.coordinates, dtype=float, copy=True)
    cells = np.array(geometry.triangles, dtype=int, copy=True)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) == 0:
        raise ValueError("invalid node array")
    if not np.all(np.isfinite(nodes)) or cells.ndim != 2 or cells.shape[1] not in (3, 6):
        raise ValueError("invalid circular element data")
    if not len(cells) or np.any(cells < 0) or np.any(cells >= len(nodes)):
        raise ValueError("invalid element node index")
    if any(len(set(map(int, c))) != len(c) for c in cells):
        raise ValueError("repeated nodes in triangular connectivity")
    want_edge_nodes = 2 if cells.shape[1] == 3 else 3
    expected = {"left", "right", "top", "bottom", "hole"}
    if set(geometry.boundary_edges) != expected:
        raise ValueError("incomplete circular boundary labels")
    edges = {}
    for label in sorted(expected):
        edge = np.array(geometry.boundary_edges[label], dtype=int, copy=True)
        if edge.ndim != 2 or edge.shape[1] != want_edge_nodes or not len(edge):
            raise ValueError(f"invalid {label} boundary edge array")
        if np.any(edge < 0) or np.any(edge >= len(nodes)):
            raise ValueError(f"out-of-range {label} boundary index")
        edges[label] = edge
    names = {name: np.unique(edge) for name, edge in edges.items()}
    for key in expected:
        supplied = np.asarray(geometry.boundary_nodes.get(key, ()), dtype=int)
        if not np.array_equal(np.sort(np.unique(supplied)), names[key]):
            raise ValueError(f"inconsistent {key} boundary nodes")
    scale = max(1., float(np.max(np.ptp(nodes, axis=0))))
    rule = list(CURVED_TRIANGLE_QUADRATURE) + [(0.,0.,0.), (1.,0.,0.), (0.,1.,0.)]
    for conn in cells:
        for r, s, _ in rule:
            det = geometric_jacobian(nodes, conn, r, s)
            if not np.isfinite(det) or det <= 1e-15 * scale**2:
                raise ValueError("nonpositive circular isoparametric Jacobian")
    return Mesh(nodes=nodes, elements=cells, edges=edges, node_sets=names)
