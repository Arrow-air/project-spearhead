import math

import numpy as np
import pytest

from spearhead.sitl import SITL_HOME, ardupilot_isa_density, sitl_home_argument
from spearhead.sitl.physics import ardupilot_pitot_airspeed


class Params:
    rho = 1.021


def indicated_by_ardupilot(json_airspeed, altitude_amsl):
    """What 4.7.1 sitl_airspeed.cpp makes of the JSON value: TAS / ISA EAS2TAS."""
    return json_airspeed/math.sqrt(1.225/ardupilot_isa_density(altitude_amsl))


def test_home_is_ankara_test_field():
    assert sitl_home_argument() == "39.858200,32.546000,1125,0"
    assert SITL_HOME[2] == 1125.0


def test_isa_density_matches_standard_atmosphere():
    assert ardupilot_isa_density(0.) == pytest.approx(1.225)
    assert ardupilot_isa_density(1000.) == pytest.approx(1.1117, abs=2e-4)


@pytest.mark.parametrize("height", [0., 80., 1000.])
def test_ardupilot_indicates_fdm_equivalent_airspeed(height):
    state = np.zeros(12)
    state[2] = -height
    state[3] = 24.
    indicated = indicated_by_ardupilot(ardupilot_pitot_airspeed(state, Params), SITL_HOME[2]+height)
    assert indicated == pytest.approx(24.*math.sqrt(1.021/1.225), rel=1e-12)


def test_pitot_cone_matches_ardupilot_rolloff():
    level = np.zeros(12); level[3] = 20.
    inside = level.copy(); inside[3], inside[5] = 20*math.cos(math.radians(15)), 20*math.sin(math.radians(15))
    outside = level.copy(); outside[3], outside[5] = 20*math.cos(math.radians(55)), 20*math.sin(math.radians(55))
    backwards = level.copy(); backwards[3] = -5.
    assert ardupilot_pitot_airspeed(inside, Params) == pytest.approx(ardupilot_pitot_airspeed(level, Params))
    assert ardupilot_pitot_airspeed(outside, Params) == pytest.approx(
        ardupilot_pitot_airspeed(level, Params)*math.cos(math.radians(35)*90/70))
    assert ardupilot_pitot_airspeed(backwards, Params) == 0.
    assert ardupilot_pitot_airspeed(np.zeros(12), Params) == 0.


def test_hover_thrust_includes_ardupilot_density_compensation():
    from spearhead.sitl import ardupilot_density_ratio
    from spearhead.sitl.quad import hover_thrust_setting
    ratio = ardupilot_density_ratio(80.)
    assert ratio == pytest.approx(ardupilot_isa_density(1205.)/1.225)
    assert hover_thrust_setting(.6663, air_density_ratio=ratio) == pytest.approx(hover_thrust_setting(.6663)*ratio)
