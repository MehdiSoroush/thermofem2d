"""Real Matplotlib/CLI postprocessing smoke tests without Streamlit dependency."""
import numpy as np
from thermofem2d.cli import main
from thermofem2d.plot import make_figures
from thermofem2d.solve import solve_model


def test_flux_and_deformation_figures(tmp_path, base_config):
    s = solve_model(base_config)
    figures = make_figures(s, tmp_path / "figures")
    names = {p.name for p in figures}
    assert {"temperature.png", "displacement.png", "von_mises.png",
            "mesh.png", "heat_flux.png", "deformed_mesh.png", "overview.png"} <= names
    assert all(p.is_file() and p.stat().st_size > 1_000 for p in figures)


def test_cli_writes_traceable_heat_flux_csv(tmp_path):
    output = tmp_path / "case"
    from pathlib import Path
    input_path = Path(__file__).resolve().parents[1] / "examples" / "steady_plate.json"
    summary = main(["--config", str(input_path), "--out", str(output)])
    assert summary["version"] == "0.1.1"
    assert summary["plate_material_area_m2"] == 3.952
    q = np.loadtxt(output / "data" / "element_heat_flux.csv", skiprows=1, delimiter=",")
    assert q.shape == (summary["elements"], 4)
    assert np.isfinite(q).all()
    assert abs(np.linalg.norm(q[:, 2:], axis=1).max() - summary["max_element_heat_flux_W_per_m2"]) < 1e-7
