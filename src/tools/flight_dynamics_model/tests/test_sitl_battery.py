import json

import numpy as np
import pytest

from spearhead.aircraft import load_aircraft_config
from spearhead.controls import Control
from spearhead.propulsion import hover_trim_control, lift_propulsion
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.battery import BatteryModel
from spearhead.sitl.physics import Physics
from spearhead.sitl.protocol import encode_sensor


def test_full_pack_open_circuit_and_resistive_sag():
    battery = BatteryModel(cells=12, capacity_ah=22.)
    assert battery.open_circuit_voltage() == pytest.approx(50.4)
    assert battery.reading(100.) == {"voltage": pytest.approx(50.4-100*.045), "current": 100.}
    with pytest.raises(ValueError):
        battery.reading(-1.)


def test_consumed_charge_lowers_voltage_and_saturates_at_empty():
    battery = BatteryModel(cells=12, capacity_ah=22.)
    battery.consume(110., 360.)          # 11 Ah: half the pack
    assert battery.state_of_charge == pytest.approx(.5)
    assert battery.open_circuit_voltage() == pytest.approx(12*3.80)
    battery.consume(1000., 3600.)
    assert battery.state_of_charge == 0 and battery.open_circuit_voltage() == pytest.approx(12*3.30)


def test_pack_current_is_lift_motors_plus_avionics_and_excludes_pusher():
    params = load_aircraft_config("pt15")
    battery = BatteryModel.from_params(params)
    hover = hover_trim_control(params)
    _, _, samples = lift_propulsion(hover, params)
    expected = 1.5 + sum(sample["current_a"] for sample in samples)
    assert battery.pack_current(hover, params) == pytest.approx(expected)
    with_pusher = Control(**{**hover.__dict__, "throttle": 1.})
    assert battery.pack_current(with_pusher, params) == pytest.approx(expected)
    assert battery.pack_current(Control(0, 0, 0, 0), params) == pytest.approx(1.5)


def test_physics_reports_battery_to_ardupilot_json_and_integrates_charge():
    params = load_aircraft_config("pt15")
    state = np.zeros(12); state[2] = -20
    battery = BatteryModel.from_params(params)
    physics = Physics(params, state, hold=False, battery=battery)
    control = hover_trim_control(params)
    sample = physics.advance(.005, control)
    current = battery.pack_current(control, params)
    assert sample["battery"]["current"] == pytest.approx(current)
    assert sample["battery"]["voltage"] == pytest.approx(battery.open_circuit_voltage()-current*.045)
    assert battery.consumed_ah == pytest.approx(current*.005/3600)
    decoded = json.loads(encode_sensor(sample))
    assert set(decoded["battery"]) == {"voltage", "current"}


def test_sitl_setup_attaches_lift_pack():
    physics, _, _ = load_setup(CONFIG_ROOT / "quad_hover.yaml")
    assert physics.battery is not None and physics.battery.cells == 12
    assert physics.battery.capacity_ah == 22.
