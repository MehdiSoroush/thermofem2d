"""Optional end-to-end Streamlit script check; executed in the UI CI job."""
from pathlib import Path
import pytest

streamlit = pytest.importorskip('streamlit')
from streamlit.testing.v1 import AppTest


def test_ui_starts_and_can_solve():
    app = Path(__file__).resolve().parents[1] / 'app.py'
    at = AppTest.from_file(str(app), default_timeout=120).run()
    assert not at.exception, 'Streamlit raised an exception during initial render'
    assert len(at.button) > 0
    at.button[0].click().run()
    assert not at.exception, 'Streamlit raised an exception while solving'
    assert len(at.success) >= 1
