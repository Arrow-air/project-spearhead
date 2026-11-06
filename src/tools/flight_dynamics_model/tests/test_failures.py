import numpy as np
import pytest
import yaml

from spearhead.controls import Control
from spearhead.propulsion import hover_trim_control
from spearhead.sitl import failures
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.failures import Failure


def hover_setup(tmp_path, failure):
    config = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    config["failure"] = failure
    path = tmp_path / "failure.yaml"; path.write_text(yaml.safe_dump(config))
    return load_setup(path)


def test_motor_out_starts_after_release_delay(tmp_path):
    physics, _, manifest = hover_setup(tmp_path, {"kind": "lift_fr_out", "start_s": .02})
    assert manifest["failure"]["kind"] == "lift_fr_out"
    hover = hover_trim_control(physics.params)
    assert physics.failure.apply(hover, physics) is hover          # still on hold
    physics.hold = False
    physics.advance(.005, hover)
    assert physics.failure.apply(hover, physics) is hover          # before start_s
    for _ in range(5):
        physics.advance(.005, hover)
    failed = physics.failure.apply(hover, physics)
    assert failed.lift_motor_throttles[0] == 0
    assert failed.lift_motor_throttles[1:] == hover.lift_motor_throttles[1:]


def test_pusher_cut_zeroes_forward_throttle(tmp_path):
    physics, _, _ = hover_setup(tmp_path, {"kind": "pusher_out", "start_s": 0})
    physics.hold = False
    physics.advance(.005, Control(0, 0, 0, .8, lift_motor_throttles=(.6,)*4))
    assert physics.failure.apply(Control(0, 0, 0, .8), physics).throttle == 0


def test_low_battery_failure_is_gated_until_parameters_are_approved(monkeypatch):
    with pytest.raises(ValueError, match="approval"):
        Failure("battery_low")
    monkeypatch.setattr(failures, "BATTERY_FAILSAFE_PARAMETERS_APPROVED", True)
    failure = Failure("battery_low", 0)

    class Pack:
        capacity_ah, consumed_ah = 22., 1.

    class FakePhysics:
        release_time, time, battery = 0., 1., Pack()

    failure.apply(Control(0, 0, 0, 0), FakePhysics)
    assert FakePhysics.battery.consumed_ah == pytest.approx(22*.85)


def test_invalid_failures_rejected():
    with pytest.raises(ValueError):
        Failure("wing_off")
    with pytest.raises(ValueError):
        Failure.from_config({"kind": "pusher_out", "when": 3})
    assert Failure.from_config(None) is None
    np.testing.assert_equal(Failure.from_config({"kind": "lift_rl_out"}).start_s, 20.)
