import subprocess
import sys

import pytest


def test_gui_app_helpers_do_not_require_nicegui():
    from spearhead.gui.app import build_parser, get_example_scenarios

    scenarios = get_example_scenarios()
    args = build_parser().parse_args(["--host", "0.0.0.0", "--port", "9090", "--no-show"])

    assert any(path.name == "pitch_doublet.yaml" for path in scenarios)
    assert args.host == "0.0.0.0"
    assert args.port == 9090
    assert args.no_show is True


def test_gui_app_registers_with_nicegui_when_available():
    pytest.importorskip("nicegui")

    from spearhead.gui.app import create_app

    create_app()


def test_gui_module_help_smoke():
    completed = subprocess.run(
        [sys.executable, "-m", "spearhead.gui", "--help"],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "python -m spearhead.gui" in completed.stdout
