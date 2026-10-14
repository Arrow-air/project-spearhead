import math

import pytest

from spearhead.gui.flight_values import VALUE_LABELS, numeric_flight_values, flight_values_status


def test_numeric_values_transform_body_velocity_and_preserve_motor_order():
    view = {'truth': {'state': [0, 0, -80, 10, 0, 0, 0, math.pi/6, -math.pi/2, .1, -.2, .3],
                      'control': {'throttle': .42}, 'pwm': [1500]*4+[1601, 1602, 1603, 1604]}}
    values = numeric_flight_values(view)
    assert values['airspeed'] == '10.00' and values['altitude'] == '80.00'
    assert values['climb'] == '5.00'  # positive upward; body w=0 is not earth-vertical speed
    assert values['heading'] == '270.0'
    view['truth']['state'][8] = math.radians(-.03)
    assert numeric_flight_values(view)['heading'] == '0.0'
    assert float(values['roll_rate']) == pytest.approx(5.7)
    assert float(values['pitch_rate']) == pytest.approx(-11.5)
    assert values['pusher'] == '42.0'
    assert [values[key] for key in ('lift_fr','lift_rr','lift_rl','lift_fl')] == ['1601','1604','1602','1603']


def test_missing_and_zero_signal_values_are_not_invented_measurements():
    assert set(numeric_flight_values({}).values()) == {'—'}
    assert list(numeric_flight_values({})) == list(VALUE_LABELS)
    row = {'state': [0]*12, 'control': {'throttle': 0}, 'pwm': [0]*8}
    values = numeric_flight_values({'truth': row})
    assert values['pusher'] == '0.0' and values['airspeed'] == '0.00'
    assert values['lift_fr'] == '—'
    row['state'][3] = float('nan')
    assert numeric_flight_values({'truth': row})['airspeed'] == '—'


def test_live_bootstrap_saved_and_stale_readings_have_distinct_status():
    assert flight_values_status({}) == ('WAITING FOR FLIGHT DATA', False)
    view = {'truth': {'hold': True}, 'active': True}
    assert flight_values_status(view) == ('STARTUP · AIRCRAFT HELD', False)
    view['truth']['hold'] = False
    assert flight_values_status(view) == ('LIVE', False)
    assert flight_values_status(dict(view, stale_truth=True))[1] is True
    assert flight_values_status(dict(view, stale_mavlink=True))[1] is True
    assert flight_values_status(dict(view, active=False)) == ('RECORDED · LAST SAMPLE', False)
