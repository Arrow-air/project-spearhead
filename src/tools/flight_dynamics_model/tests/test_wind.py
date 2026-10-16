import math

import numpy as np
import pytest
import yaml

from spearhead.aircraft import load_aircraft_config
from spearhead.controls import Control
from spearhead.force_moment import compute_force_moment_breakdown
from spearhead.propulsion import hover_trim_control
from spearhead.sitl.backend import CONFIG_ROOT, load_setup
from spearhead.sitl.physics import Physics, ardupilot_pitot_airspeed
from spearhead.wind import Gust, WindField


def test_meteorological_direction_and_one_cosine_gust():
    np.testing.assert_allclose(WindField(5., 0.).ned(3.), [-5, 0, 0], atol=1e-12)     # from north
    np.testing.assert_allclose(WindField(5., 90.).ned(0.), [0, -5, 0], atol=1e-12)    # from east
    gust = Gust(start_s=10., duration_s=4., speed_m_s=3., from_deg=0.)
    field = WindField(0., 0., (gust,))
    assert np.allclose(field.ned(9.9), 0) and np.allclose(field.ned(14.1), 0)
    np.testing.assert_allclose(field.ned(12.), [-3, 0, 0], atol=1e-12)
    np.testing.assert_allclose(field.ned(11.), [-1.5, 0, 0], atol=1e-12)
    with pytest.raises(ValueError):
        WindField(25.)
    with pytest.raises(ValueError):
        Gust(0., 0., 1., 0.)


def test_config_round_trip_and_validation():
    config = {"speed_m_s": 6, "from_deg": 270, "gusts": [{"start_s": 20, "duration_s": 3, "speed_m_s": 4, "from_deg": 270}]}
    field = WindField.from_config(config)
    assert WindField.from_config(field.as_config()) == field
    assert WindField.from_config(None) is None and WindField().calm
    with pytest.raises(ValueError):
        WindField.from_config({"speed": 5})


def test_aircraft_moving_with_the_air_has_no_aero_or_axial_inflow():
    params = load_aircraft_config("pt15")
    x = np.zeros(12); x[3] = 8.
    calm = compute_force_moment_breakdown(0, x, params, Control(0, 0, 0, 0))
    # A tailwind equal to the ground speed leaves zero airspeed.
    drifting = compute_force_moment_breakdown(0, x, params, Control(0, 0, 0, 0), wind_ned=np.array([8., 0, 0]))
    still = compute_force_moment_breakdown(0, np.zeros(12), params, Control(0, 0, 0, 0))
    np.testing.assert_allclose(drifting["wing"]["force_b"], still["wing"]["force_b"], atol=1e-9)
    assert np.linalg.norm(calm["wing"]["force_b"]) > 1
    # A headwind over a hovering aircraft produces the same aero as flying into still air.
    headwind = compute_force_moment_breakdown(0, np.zeros(12), params, Control(0, 0, 0, 0), wind_ned=np.array([-8., 0, 0]))
    np.testing.assert_allclose(headwind["wing"]["force_b"], calm["wing"]["force_b"], rtol=1e-9)


def test_pitot_reads_air_relative_speed():
    params = load_aircraft_config("pt15")
    x = np.zeros(12)
    still = ardupilot_pitot_airspeed(np.r_[0, 0, 0, 10, 0, 0, np.zeros(6)], params)
    assert ardupilot_pitot_airspeed(x, params, np.array([-10., 0, 0])) == pytest.approx(still)
    assert ardupilot_pitot_airspeed(x, params, np.array([10., 0, 0])) == 0     # tailwind: pitot sees nothing


def test_physics_holds_still_air_until_release_then_applies_wind():
    params = load_aircraft_config("pt15")
    state = np.zeros(12); state[2] = -20
    wind = WindField(6., 0.)
    physics = Physics(params, state, hold=True, wind=wind)
    hover = hover_trim_control(params)
    assert physics.advance(.005, hover)["airspeed"] == 0
    assert physics.wind_ned(physics.time) is None
    physics.hold = False
    sample = physics.advance(.005, hover)
    np.testing.assert_allclose(physics.wind_ned(physics.time), [-6, 0, 0])
    assert sample["airspeed"] > 5
    # The headwind's drag acts backwards on the free aircraft.
    loads = compute_force_moment_breakdown(physics.time, physics.state, params, hover,
                                           wind_ned=physics.wind_ned(physics.time))
    assert (loads["wing"]["force_b"] + loads["fuselage"]["force_b"] + loads["tail"]["force_b"])[0] < 0


def test_sitl_config_accepts_wind_and_records_it(tmp_path):
    config = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    config["wind"] = {"speed_m_s": 5, "from_deg": 45, "gusts": [{"start_s": 15, "duration_s": 3, "speed_m_s": 3, "from_deg": 45}]}
    path = tmp_path / "windy.yaml"; path.write_text(yaml.safe_dump(config))
    physics, _, manifest = load_setup(path)
    assert physics.wind.speed_m_s == 5 and len(physics.wind.gusts) == 1
    assert manifest["wind"]["from_deg"] == 45
    calm = yaml.safe_load((CONFIG_ROOT / "quad_hover.yaml").read_text())
    path.write_text(yaml.safe_dump(calm))
    physics, _, manifest = load_setup(path)
    assert physics.wind is None and manifest["wind"] is None
