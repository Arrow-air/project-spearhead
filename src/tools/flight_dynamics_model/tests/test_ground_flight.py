from dataclasses import replace

import numpy as np
import pytest

from spearhead.controls import Control
from spearhead.propulsion import hover_trim_control
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.physics import sensor_state


def test_ground_equilibrium_and_accelerometer():
    physics, _, _ = load_setup(CONFIG_ROOT / "quad_ground.yaml")
    off = Control(0, 0, 0, 0)
    force, moment, contact = physics.ground.loads(physics.state)
    np.testing.assert_allclose(force, [0, 0, -physics.params.mass*physics.params.g], atol=1e-10)
    np.testing.assert_allclose(moment, 0, atol=1e-10)
    assert contact["contact"] and not physics.hold
    np.testing.assert_allclose(physics.derivative(0, physics.state, off), 0, atol=1e-10)
    sensor = sensor_state(0, physics.state, physics.params, off, ground=physics.ground)
    np.testing.assert_allclose(sensor["imu"]["accel_body"], [0, 0, -physics.params.g], atol=1e-10)


def test_liftoff_and_landing_are_free_physics():
    physics, _, _ = load_setup(CONFIG_ROOT / "quad_ground.yaml")
    rest = physics.state.copy()
    control = hover_trim_control(replace(physics.params, mass=physics.params.mass * 1.2))
    for _ in range(200):
        physics.advance(.005, control)
    assert -physics.state[2] > 1
    assert not physics.ground.loads(physics.state)[2]["contact"]
    # Release from 0.6 m with motors off; the gear must absorb the drop.
    physics.state = rest.copy(); physics.state[2] = -.6
    for _ in range(400):
        physics.advance(.005, Control(0, 0, 0, 0))
    assert physics.ground.loads(physics.state)[2]["contact"]
    np.testing.assert_allclose(physics.state, rest, atol=.001)


def test_no_ground_loads_in_air_and_friction_is_bounded():
    physics, _, _ = load_setup(CONFIG_ROOT / "quad_ground.yaml")
    state = physics.state.copy(); state[2] = -10
    force, moment, contact = physics.ground.loads(state)
    np.testing.assert_array_equal(force, 0); np.testing.assert_array_equal(moment, 0)
    assert not contact["contact"]
    state = physics.state.copy(); state[3] = 100
    force, _, contact = physics.ground.loads(state)
    assert force[0] < 0
    assert abs(force[0]) <= physics.ground.friction*contact["normal_force_n"]+1e-9


def test_subweight_lift_stays_grounded():
    physics, _, manifest = load_setup(CONFIG_ROOT / "quad_ground.yaml")
    control = hover_trim_control(replace(physics.params, mass=physics.params.mass * .8))
    for _ in range(100):
        physics.advance(.005, control)
    assert physics.ground.loads(physics.state)[2]["contact"]
    assert abs(physics.state[5]) < .001


def test_air_start_still_stops_at_ground():
    physics, _, _ = load_setup(CONFIG_ROOT / "quad_hover.yaml")
    physics.hold = False; physics.state[2] = -.001; physics.state[5] = 1
    with pytest.raises(ValueError):
        physics.advance(.005, Control(0, 0, 0, 0))


def test_mission_upload_requires_every_item_and_acceptance():
    from types import SimpleNamespace
    from spearhead.sitl.ground_flight import GroundMission
    calls = []
    mav = SimpleNamespace(mission_count_send=lambda *a: calls.append(a),
                          mission_item_int_send=lambda *a: calls.append(a))
    link = SimpleNamespace(mav=mav, target_system=1, target_component=1)
    mission = GroundMission(5, 15)
    mission.update(link, None, 0)
    assert calls[0][-1] == 5
    ack = SimpleNamespace(type=0, get_type=lambda: "MISSION_ACK")
    with pytest.raises(RuntimeError):
        mission.update(link, ack, 1)
    for seq in range(5):
        request = SimpleNamespace(seq=seq, get_type=lambda: "MISSION_REQUEST")
        mission.update(link, request, 2)
    assert not mission.accepted
    mission.update(link, ack, 3)
    assert mission.accepted
    # Explicit MC command is essential: loiter must never initiate a FW transition.
    assert calls[3][4] == 3000 and calls[3][7] == 3
