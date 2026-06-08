import numpy as np

from spearhead.force_moment import COMPONENTS, compute_force_moment_breakdown
from spearhead.params import AircraftParams
from spearhead.simulation import default_initial_state


def test_force_moment_breakdown_has_expected_components_and_shapes():
    params = AircraftParams()
    breakdown = compute_force_moment_breakdown(0.0, default_initial_state(), params)

    assert set(breakdown) == set(COMPONENTS)
    for component in COMPONENTS:
        assert breakdown[component]["force_b"].shape == (3,)
        assert breakdown[component]["moment_b"].shape == (3,)
        assert np.all(np.isfinite(breakdown[component]["force_b"]))
        assert np.all(np.isfinite(breakdown[component]["moment_b"]))


def test_force_moment_breakdown_total_matches_component_sum():
    params = AircraftParams()
    breakdown = compute_force_moment_breakdown(5.0, default_initial_state(), params)
    non_total = [component for component in COMPONENTS if component != "total"]

    force_sum = sum((breakdown[component]["force_b"] for component in non_total), np.zeros(3))
    moment_sum = sum((breakdown[component]["moment_b"] for component in non_total), np.zeros(3))

    np.testing.assert_allclose(breakdown["total"]["force_b"], force_sum)
    np.testing.assert_allclose(breakdown["total"]["moment_b"], moment_sum)


def test_aero_component_split_sums_to_total_aero_load():
    params = AircraftParams()
    breakdown = compute_force_moment_breakdown(11.0, default_initial_state(), params)
    aero_components = ("wing", "tail", "fuselage")

    aero_force_sum = sum((breakdown[component]["force_b"] for component in aero_components), np.zeros(3))
    aero_moment_sum = sum((breakdown[component]["moment_b"] for component in aero_components), np.zeros(3))
    non_aero_force = breakdown["propulsion"]["force_b"] + breakdown["gravity"]["force_b"]
    non_aero_moment = breakdown["propulsion"]["moment_b"] + breakdown["gravity"]["moment_b"]

    np.testing.assert_allclose(aero_force_sum, breakdown["total"]["force_b"] - non_aero_force)
    np.testing.assert_allclose(aero_moment_sum, breakdown["total"]["moment_b"] - non_aero_moment)

