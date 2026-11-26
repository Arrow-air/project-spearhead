"""PT1.5 lift curve: MAD V8013 150KV + AMPX 80A + FLUXER PRO 26x7.8 at Ankara density."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from spearhead import load_aircraft_config
from spearhead.cg_card import hover_headroom
from spearhead.controls import Control
from spearhead.propulsion import hover_collective, lift_propulsion, load_static_performance

FDM_ROOT = Path(__file__).resolve().parents[1]
PROPULSION = FDM_ROOT / "data/propulsion"


def _builder():
    spec = importlib.util.spec_from_file_location("build_mad_lift_curve", FDM_ROOT / "scripts/build_mad_lift_curve.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_source_matches_the_stored_sheet_and_is_self_consistent():
    source = json.loads((PROPULSION / "mad_v8013_150_fluxer26_source.json").read_text())
    image = PROPULSION / source["source_file"]
    assert hashlib.sha256(image.read_bytes()).hexdigest() == source["source_file_sha256"]
    table = np.asarray(source["samples"])
    throttle, voltage, current, power_in, power_out, torque, rpm, thrust, efficiency, gf_per_w = table.T
    assert len(table) == 15 and throttle[0] == .3 and throttle[-1] == 1
    np.testing.assert_allclose(power_in, voltage*current, rtol=.002)
    np.testing.assert_allclose(power_out, torque*rpm*2*np.pi/60, rtol=.003)
    np.testing.assert_allclose(efficiency, 100*power_out/power_in, atol=.11)
    np.testing.assert_allclose(gf_per_w, thrust/power_in, atol=.051)
    # The older pusher-fit transcription carries the same rows without thrust; it reads
    # 3810 RPM at 60 % where the sheet shows 3811 (and output power confirms 3811).
    bench = json.loads((PROPULSION / "mad_v8013_150_bench.json").read_text())
    np.testing.assert_allclose(np.asarray(bench["samples"]), table[:, [0, 1, 2, 5, 6]], rtol=0, atol=1)
    assert power_out[6]/(torque[6]*2*np.pi/60) == pytest.approx(3811, abs=.5)


def test_rebuilt_curve_reproduces_the_sheet_and_the_stored_files(tmp_path):
    builder = _builder()
    data, sheet = builder.load_sheet()
    model = builder.fit(data, sheet)
    _, summary = builder.validation(model, sheet)
    assert summary["thrust_error"] < .06 and summary["current_error"] < .08 and summary["rpm_error"] < .03
    stored = json.loads((PROPULSION / "mad_v8013_150_fluxer26_fit.json").read_text())
    assert stored["equivalent_resistance_ohm"] == pytest.approx(model["resistance"], rel=1e-6)
    assert stored["curve_sha256"] == hashlib.sha256((PROPULSION / stored["curve"]).read_bytes()).hexdigest()
    assert stored["source_sha256"] == hashlib.sha256(builder.SOURCE.read_bytes()).hexdigest()
    # Density scaling: same throttle and voltage at 1.021 kg/m3 gives less thrust, slightly more RPM.
    sea = builder.operating_point(model, 48., .65, model["test_density"])
    ankara = builder.operating_point(model, 48., .65, 1.021)
    assert ankara["thrust_n"] < .9*sea["thrust_n"] and ankara["rpm"] > sea["rpm"]
    assert ankara["current_a"] < sea["current_a"]


def test_profile_uses_the_curve_over_the_12s_range():
    params = load_aircraft_config("pt15")
    assert params.lift_system.performance_path.name == "mad_v8013_150_fluxer26_rho1021.csv"
    curve = load_static_performance(params.lift_system.performance_path)
    assert curve.voltages[0] == 36 and curve.voltages[-1] == 50.4
    assert curve.throttles[1] == .15 and curve.throttles[-1] == 1
    with (PROPULSION / "mad_v8013_150_fluxer26_rho1021.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    row = next(r for r in rows if r["voltage_v"] == "44.4" and r["throttle"] == "0.65")
    point = curve.lookup(44.4, .65)
    assert point["thrust_n"] == pytest.approx(float(row["thrust_n"]))
    assert {r["basis"] for r in rows if float(r["throttle"]) < .3} == {"model_below_sheet"}
    # Full throttle at 44.4 V: about 12.8 kgf per rotor at Ankara (15.0 kgf on the 48 V sea-level sheet).
    assert 120 < curve.lookup(44.4, 1)["thrust_n"] < 130


def test_hover_at_22_kg_stays_inside_the_sustainable_motor_current():
    params = load_aircraft_config("pt15")
    collective = hover_collective(params)
    assert .62 < collective < .67
    _, _, samples = lift_propulsion(Control(0, 0, 0, 0, lift_throttle=collective), params)
    assert sum(s["thrust_n"] for s in samples) == pytest.approx(params.mass*params.g)
    current = [s["current_a"] for s in samples]
    assert all(17 < value < params.lift_system.motor["sustainable_current_a"] for value in current)
    assert all(s["torque_nm"] > 1.5 for s in samples)
    # Busiest rotor against its 95 % spin thrust at 44.4 V (KDE PT1.5 at 21 kg: 65 %).
    assert .45 < hover_headroom(params) < .5
