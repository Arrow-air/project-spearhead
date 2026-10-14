import json

import pytest

from spearhead.gui.flight_setup import FlightSetup, AIRCRAFT
from spearhead.sitl.transition import analyse_transition, transition_command


def test_forward_trim_preserves_hover_rc_mapping_and_actual_prototype_cg():
    from spearhead.sitl.backend import load_setup, CONFIG_ROOT
    from spearhead.sitl.transition import forward_trim_parameters
    physics, mapping, manifest = load_setup(CONFIG_ROOT/'quad_transition.yaml')
    cg = physics.params.cg_body_xyz.copy()
    parameters, _, control = forward_trim_parameters(physics.params, mapping, 80)
    assert list(physics.params.cg_body_xyz) == list(cg) == [.01, 0, .15]
    assert 'RC3_TRIM' not in parameters and 'Q_M_THST_HOVER' not in parameters
    # Reference correction changes the solved elevator sign, not its mapping.
    pwm = [1500] * 16
    pwm[mapping.elevator.number - 1] = parameters['SERVO2_TRIM']
    assert abs(mapping.control(pwm).de - control.de) < .00053
    assert isinstance(parameters['TRIM_THROTTLE'], int)
    assert abs(control.de) < .5236 and control.lift_throttle == 0


def test_forward_throttle_waits_for_confirmed_fbwa():
    assert transition_command(9, 18) == (18, 1500)
    assert transition_command(15, 18) == (5, 1500)
    assert transition_command(20, 5) == (5, 1500)  # No confirmation timestamp.
    assert transition_command(12.5, 5, 10) == (5, 1700)
    assert transition_command(20, 5, 10) == (5, 1900)


def test_delayed_confirmation_gets_the_full_five_second_ramp():
    assert transition_command(15, 18) == (5, 1500)
    assert transition_command(15.01, 5, 15.01) == (5, 1500)
    assert transition_command(17.51, 5, 15.01) == (5, 1700)
    assert transition_command(20.01, 5, 15.01) == (5, 1900)
    assert transition_command(18, 18, 15.01) == (5, 1500)
    assert transition_command(19, 5, 19) == (5, 1500)


def test_transition_requires_prototype_air_start_and_altitude_margin():
    assert AIRCRAFT['pt1_as_built'] == 'PT1.5 · 21 kg'
    setup = FlightSetup(initial_point='hover80', maneuver='transition', altitude_m=80)
    assert setup.validate().scenario == 'quad_transition'
    assert setup.config()['cg_body_xyz_m'] == [.01, 0, .15]
    with pytest.raises(ValueError, match='80 m'):
        FlightSetup(maneuver='transition').validate()
    with pytest.raises(ValueError, match='maneuver'):
        FlightSetup(aircraft='reference_design', initial_point='trim28', maneuver='transition', speed_m_s=28, altitude_m=100).validate()


def test_mode_change_and_speed_are_insufficient_if_lift_motors_still_support_flight(tmp_path):
    path = tmp_path/'truth.jsonl'
    def write(lift_pwm):
        rows = [{'time_s': t, 'hold': False, 'state': [0,0,-80,25,0,0,0,0,0,0,0,0],
                 'pwm': [1500,1500,1900,1500]+[lift_pwm]*4} for t in range(61)]
        path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        (tmp_path/'telemetry.jsonl').write_text(''.join(json.dumps({
            'fdm_time_s':t, 'message':{'mavpackettype':'HEARTBEAT','custom_mode':5}})+'\n'
            for t in range(61)))
    write(1740)
    report = analyse_transition(path, 0, 60, 80, 5)
    assert not report['passed'] and report['longest_forward_without_lift_s'] == 0
    write(1100)
    assert analyse_transition(path, 0, 60, 80, 5)['passed']
    write(1185)
    assert analyse_transition(path, 0, 60, 80, 5)['passed']
    write(1250)
    assert not analyse_transition(path, 0, 60, 80, 5)['passed']
    assert not analyse_transition(path, 0, 60, 80, 18)['passed']


def test_final_fbwa_does_not_confirm_an_earlier_interval_or_missing_heartbeats(tmp_path):
    path = tmp_path/'truth.jsonl'
    rows = [{'time_s':t, 'hold':False, 'state':[0,0,-80,25,0,0,0,0,0,0,0,0],
             'pwm':[1500,1500,1900,1500]+[1100]*4} for t in range(61)]
    path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
    assert not analyse_transition(path,0,60,80,5)['passed']
    telemetry = tmp_path/'telemetry.jsonl'
    def modes(values):
        telemetry.write_text(''.join(json.dumps({'fdm_time_s':t,
            'message':{'mavpackettype':'HEARTBEAT','custom_mode':mode}})+'\n'
            for t,mode in values))
    modes([(0,5),(60,5)])  # Stale confirmation cannot span a missing log.
    assert not analyse_transition(path,0,60,80,5)['passed']
    modes([(t,5 if t >= 55 else 18) for t in range(61)])
    result = analyse_transition(path,0,60,80,5)
    assert not result['passed'] and result['longest_forward_without_lift_s'] == 5
    modes([(t,5 if t % 8 != 0 else 18) for t in range(61)])
    assert not analyse_transition(path,0,60,80,5)['passed']
    modes([(t,5) for t in range(61)])
    assert analyse_transition(path,0,60,80,5)['passed']


def test_forward_mode_join_never_uses_future_confirmation(tmp_path):
    import numpy as np
    from spearhead.sitl.transition import confirmed_forward_samples
    path = tmp_path/'telemetry.jsonl'
    path.write_text(json.dumps({'fdm_time_s':10,
        'message':{'mavpackettype':'HEARTBEAT','custom_mode':5}})+'\n')
    mask, evidence = confirmed_forward_samples(np.array([9.99,10.,12.,12.01]),path)
    assert evidence and mask.tolist() == [False,True,True,False]
