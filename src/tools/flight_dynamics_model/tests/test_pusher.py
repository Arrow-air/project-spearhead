"""Check units, loaded motor equilibrium and state-dependent FDM integration."""

from dataclasses import replace
import hashlib
import json

import numpy as np
import pytest

from spearhead import load_aircraft_config
from spearhead.aircraft import AIRCRAFT_CONFIGS
from spearhead.controls import Control
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.propulsion import pusher_performance, pusher_propulsion, propulsion_model
from spearhead.pusher import estimate_pusher, load_propeller_map


@pytest.fixture
def pt1():
    return load_aircraft_config("pt1_as_built")


def test_motor_fit_matches_manufacturer_operating_points(pt1):
    source = json.loads((AIRCRAFT_CONFIGS.parent / "data/propulsion/mad_v8013_150_bench.json").read_text())
    samples = np.asarray(source["samples"])
    duty, voltage, _, torque, rpm = samples.T
    kt = 60/(2*np.pi*150)
    winding_current = torque/kt+2.35
    required_voltage = rpm/150+pt1.pusher_system.equivalent_resistance_ohm*winding_current
    # Fit checked against all 15 points, not the much larger 26-inch prop thrust.
    assert np.sqrt(np.mean((required_voltage-duty*voltage)**2)) < .21
    assert np.max(abs(required_voltage-duty*voltage)) < .6
    metadata = json.loads((AIRCRAFT_CONFIGS.parent / "data/propulsion/pt1_pusher_source.json").read_text())
    for source in metadata["files"]:
        path = AIRCRAFT_CONFIGS.parent / "data/propulsion" / source["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source["sha256"]


def test_apc_coefficients_reconstruct_published_dimensional_units(pt1):
    propeller = load_propeller_map(pt1.pusher_system.performance_path)
    # APC's 1000 RPM static row: Ct=.0838, Cp=.0337, 1.900 N, 6.468 W.
    ct, cp, outside = propeller.coefficients(1000, 0)
    assert not outside
    diameter = 20*.0254
    n = 1000/60
    assert ct*1.225*n**2*diameter**4 == pytest.approx(1.900, rel=.002)
    assert cp*1.225*n**3*diameter**5 == pytest.approx(6.468, rel=.002)
    assert propeller.coefficients(1500, 0)[0] == pytest.approx((.0838+.0842)/2)


@pytest.mark.parametrize("speed", [0, 10, 20, 22, 24, 26, 28])
def test_estimate_respects_voltage_torque_and_energy_balance(pt1, speed):
    system = pt1.pusher_system
    sample = estimate_pusher(system, 44.4, .85, speed, 1.021)
    kt = 60/(2*np.pi*150)
    winding_current = sample["torque_nm"]/kt+2.35
    assert sample["rpm"]/150+winding_current*system.equivalent_resistance_ohm == pytest.approx(.85*44.4)
    assert sample["power_w"] == pytest.approx(sample["current_a"]*44.4)
    assert sample["shaft_power_w"] == pytest.approx(sample["torque_nm"]*sample["rpm"]*2*np.pi/60)
    assert 0 < sample["shaft_power_w"] < sample["power_w"]
    assert sample["rpm"] < .85*44.4*150
    assert not sample["outside_propeller_data"]
    json.dumps(sample)  # Telemetry needs native JSON types.


def test_maximum_limits_and_sustainable_guideline_are_distinct(pt1):
    sample = estimate_pusher(pt1.pusher_system, 44.4, 1, 0, 1.021)
    assert sample["above_sustainable_current"] and not sample["limited"]
    constrained = replace(pt1.pusher_system, max_current_a=24)
    sample = estimate_pusher(constrained, 44.4, 1, 10, 1.021)
    assert sample["current_a"] == pytest.approx(24)
    assert sample["limited"] and sample["effective_duty"] < 1
    constrained = replace(pt1.pusher_system, max_input_power_w=600)
    sample = estimate_pusher(constrained, 44.4, 1, 0, 1.021)
    assert sample["power_w"] == pytest.approx(600)


def test_static_and_forward_estimates_use_correct_state_and_hardware(pt1):
    control = Control(0, 0, 0, 1)
    static = pusher_performance(control, pt1)
    forward = pusher_performance(control, pt1, 26)
    json.dumps(pusher_performance(control, pt1, np.float64(26)))
    assert static["thrust_n"] == pytest.approx(pt1.Tmax)
    assert forward["thrust_n"] < static["thrust_n"]
    state = np.zeros(12)
    state[3] = 26
    loads = compute_force_moment_breakdown(0, state, pt1, control)
    assert loads["propulsion"]["force_b"][0] == pytest.approx(forward["thrust_n"])
    reference = load_aircraft_config("reference_design")
    assert propulsion_model(control, reference, 26)[0][0] == 120  # Reference unchanged.
    assert pusher_performance(control, replace(pt1, pusher_installed=False))["thrust_n"] == 0


def test_off_bounds_and_unmodelled_inflow_are_explicit(pt1):
    system = pt1.pusher_system
    assert estimate_pusher(system, 44.4, 0, 26, 1.021)["power_w"] == 0
    assert estimate_pusher(system, 44.4, -1, 0, 1.021)["thrust_n"] == 0
    assert estimate_pusher(system, 44.4, 2, 0, 1.021) == estimate_pusher(system, 44.4, 1, 0, 1.021)
    for args in [(60, .5, 0, 1.021), (44.4, float("nan"), 0, 1.021), (44.4, .5, 0, -1)]:
        with pytest.raises(ValueError):
            estimate_pusher(system, *args)
    assert estimate_pusher(system, 44.4, .01, 0, 1.021)["outside_propeller_data"]
    assert estimate_pusher(system, 44.4, .5, -5, 1.021)["outside_propeller_data"]
    high_advance = estimate_pusher(system, 44.4, .5, 100, 1.021)
    assert high_advance["outside_propeller_data"] and high_advance["thrust_n"] == 0


def test_optional_mount_wrench_and_independent_lift_command(pt1):
    position = tuple(pt1.cg_body_xyz + np.array([-.8, .1, -.05]))
    pt1 = replace(pt1, pusher_system=replace(pt1.pusher_system, position_body_m=position, reaction_torque_sign=-1))
    control = Control(0, 0, 0, .6)
    force, moment = pusher_propulsion(control, pt1, 22)
    sample = pusher_performance(control, pt1, 22)
    expected = np.cross(np.asarray(position)-pt1.cg_body_xyz, force)
    expected[0] -= sample["torque_nm"]
    np.testing.assert_allclose(moment, expected)
    combined, _ = propulsion_model(replace(control, lift_throttle=.5), pt1, 22)
    assert combined[0] == force[0] and combined[2] < 0
