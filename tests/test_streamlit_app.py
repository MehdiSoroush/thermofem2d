"""Optional end-to-end Streamlit script check; executed in the UI CI job."""
from pathlib import Path
import pytest

streamlit = pytest.importorskip('streamlit')
from streamlit.testing.v1 import AppTest


def test_ui_starts_and_can_solve():
    app = Path(__file__).resolve().parents[1] / 'app.py'
    at = AppTest.from_file(str(app), default_timeout=120).run()
    assert not at.exception, 'Streamlit raised an exception during initial render'
    # Button is in the sidebar. Run a real numerical model through the UI.
    assert len(at.button) > 0
    at.button[0].click().run()
    assert not at.exception, 'Streamlit raised an exception while solving'
    assert len(at.success) >= 1


def test_ui_can_solve_with_t6():
    app=Path(__file__).resolve().parents[1]/'app.py'
    at=AppTest.from_file(str(app),default_timeout=120).run()
    assert not at.exception
    # The first selectbox in the sidebar chooses the element interpolation.
    labels=[box.label for box in at.selectbox]
    idx=labels.index('Finite element')
    at.selectbox[idx].set_value('T6').run()
    at.button[0].click().run()
    assert not at.exception
    assert any('T6' in item.message for item in at.success)
