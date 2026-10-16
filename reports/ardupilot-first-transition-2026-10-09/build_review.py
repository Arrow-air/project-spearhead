"""Build a read-only proposal from logged values and pinned 4.7.1 evidence."""
import csv, hashlib, json, math, pathlib, re, shutil

root = pathlib.Path(__file__).parent
repo = root.parent.parent
ap = pathlib.Path('/private/tmp/spearhead-ardupilot-4.7.1')
commit = 'dbe792162d06cab66c3475fd5556bf7a120f119e'
logged = json.loads((root/'log17-last-parameters.json').read_text())
sitl = json.loads((root/'sitl-last-parameters.json').read_text())
manifest = json.loads((repo/'artifacts/fdm_sitl_envelope_2026_10_09/roll_plus_28/manifest.json').read_text())

def parse(text):
    result = {}
    for line in text.splitlines():
        fields = line.split('#', 1)[0].replace(',', ' ').split()
        if len(fields) == 2:
            result[fields[0]] = float(fields[1])
    return result

peers, identities = {}, {}
for name in ['Foxtech_GreatShark', 'Mugin_EV350']:
    path = ap/'Tools/Frame_params/QuadPlanes'/f'{name}.param'
    content = path.read_bytes()
    (root/f'{name}.source.param').write_bytes(content)
    peers[name] = parse(content.decode())
    identities[name] = {'source_commit':commit, 'path':str(path), 'sha256':hashlib.sha256(content).hexdigest(), 'url':f'https://github.com/ArduPilot/ardupilot/blob/{commit}/Tools/Frame_params/QuadPlanes/{name}.param', 'license':'ArduPilot GPL-3.0; upstream file retained with comments'}

sources = {
 'transitions':'https://ardupilot.org/plane/docs/quadplane-transitions.html',
 'reliability':'https://ardupilot.org/plane/docs/quadplane-reliability.html',
 'airspeed':'https://ardupilot.org/plane/docs/airspeed-parameters-setup.html',
 'sensor':'https://ardupilot.org/plane/docs/calibrating-an-airspeed-sensor.html',
 'tuning':'https://ardupilot.org/plane/docs/new-roll-and-pitch-tuning.html',
 'vtol':'https://ardupilot.org/plane/docs/quadplane-vtol-tuning-process.html',
 'tecs':'https://ardupilot.org/plane/docs/tecs-total-energy-control-system-for-speed-height-tuning-guide.html',
 'flying':'https://ardupilot.org/plane/docs/quadplane-flying.html',
 'qsource':f'https://github.com/ArduPilot/ardupilot/blob/{commit}/ArduPlane/quadplane.cpp',
 'psource':f'https://github.com/ArduPilot/ardupilot/blob/{commit}/ArduPlane/Parameters.cpp',
 'bsource':f'https://github.com/ArduPilot/ardupilot/blob/{commit}/libraries/AP_BattMonitor/AP_BattMonitor_Params.cpp',
 'greatshark':'https://www.foxtechfpv.com/foxtech-great-shark-330-vtol.html',
 'mugin':'https://www.muginuav.com/product/mugin-ev350-carbon-fiber-full-electric-vtol-uav-platform/',
 'kde':'https://www.kdedirect.com/products/kde5215xf-220',
 'esc':'https://www.kdedirect.com/products/kde-uas75hvc',
}

# No values are written to the controller, FDM runner or GUI by this script.
proposal = {}
def add(name, value, reason, evidence, condition='Bench/SITL review before aircraft loading'):
    proposal[name] = {'value':value, 'reason':reason, 'source':sources[evidence], 'condition':condition}

speed_gate = 'Illustrative 25 kg winged configuration ONLY: calibrated indicated stall speed <=14.4 m/s; replace from measured AUW/stall margin; corrected controls and real CG still required'
add('AIRSPEED_MIN',18,'Replace wingless 9 m/s; 1.25 times a 14.4 m/s indicated stall speed. Not derived from a measured stall test.','reliability',speed_gate)
add('Q_ASSIST_SPEED',15,'Enable assistance; MIN minus 3 m/s, above the assumed 14.4 m/s stall. In 4.7.1 -1 disables ALL Q_ASSIST triggers outside transitions.','qsource',speed_gate)
add('AIRSPEED_CRUISE',25,'First steady-flight target near the existing 28 m/s true-speed test after density conversion; 25 indicated is ~27.38 true at rho=1.021.','psource',speed_gate)
add('AIRSPEED_MAX',27,'Automatic demanded-speed ceiling; 1.5 times the illustrative MIN=18. It is not an overspeed limiter and is not a requested test speed.','psource',speed_gate+'; confirm attainable sustainable full-throttle level speed; no tests at 30 m/s or above')
add('SCALING_SPEED',25,'Set the speed-scaling reference near cruise before the first fixed-wing tune; logged fixed-wing gains match the 4.7.1 defaults. Existing SITL tune uses 28 and cannot be transferred.','airspeed','Only with default/reset fixed-wing rate gains and a new tune; verify physical surface mapping first')
add('Q_TRANS_FAIL',20,'Bound airspeed-wait time instead of indefinite transition; peer GreatShark also uses 20 s. Demonstrate acceleration to MIN within this time.','qsource','Validate real forward-thrust curve and timeout/abort in QuadPlane SITL')
add('Q_TRANS_FAIL_ACT',1,'QRTL on failed transition; 4.7.1 enum 1=QRTL, 0=QLAND.','qsource','Validate QRTL approach, return path, home and battery reserve; QRTL can command fixed-wing flight when distant')
add('Q_OPTIONS',int(logged['Q_OPTIONS']) | 1 | 32,'Add bit 0 level transition and bit 5 VTOL RC-failsafe QRTL; preserve bit 18 VTOL-only arming. Do not add bit 19 or 20.','qsource','Review bitmask against fresh board export and rehearse QRTL; level transition restricts bank and VTOL climb response')
add('ROLL_LIMIT_DEG',25,'Reduce maximum fixed-wing bank from 45 degrees for initial gentle manoeuvres; level-turn stall multiplier decreases from 1.189 to 1.050.','psource','First-flight operating limit, not an aerodynamic structural limit')
add('PTCH_LIM_MAX_DEG',10,'Reduce nose-up demand from 20 degrees while climb capability is unmeasured; matches present SITL angle envelope only.','tecs','Confirm sustainable climb and surface authority; FBWA does not regulate airspeed or altitude')
add('PTCH_LIM_MIN_DEG',-10,'Reduce nose-down demand from -25 degrees for the initial flight; matches SITL envelope.','tecs','Pilot retains throttle/speed responsibility in FBWA')
add('PTCH2SRV_TCONST',0.75,'Use a gentler initial pitch demand than 0.5 logged / 0.2 SITL; both similar-aircraft reference files use 0.75. This is a demand-shaping choice, not a validated tune.','tuning','Recheck actual pitch response before TECS/autonomous flight')
add('RLL2SRV_RMAX',45,'Replace unlimited roll rate demand with 45 deg/s for the initial test.','psource','Engineering starting limit; adjust from actual response')
add('PTCH2SRV_RMAX_UP',22,'Limit upward pitch rate; approximately 560/25=22.4 deg/s per tuning guidance.','tuning','Rate-demand limit only; does not certify load factor or prevent stall')
add('PTCH2SRV_RMAX_DN',22,'Apply the same initial pitch demand limit downward.','tuning','Rate-demand limit only; verify elevator direction and measured throw')
add('FLTMODE3',19,'Expose QLOITER immediately before the FBWA slot.','flying','Verify actual six-position RC7 switch mapping and pilot bailout')
add('FLTMODE4',5,'Expose FBWA for manual-forward-throttle first transition; do not start in CRUISE/AUTO.','flying','Wings and forward propulsion installed; surface/pusher direction checked')
add('FLTMODE5',21,'Expose QRTL return mode after FBWA.','flying','Rehearse the return path and confirm switch positions')
add('FLTMODE6',20,'Expose QLAND at the final switch slot.','flying','Verify mode selection on the actual radio')
add('FS_LONG_ACTN',1,'Explicit RTL on long failsafe rather than continuing an automatic mode.','psource','End-to-end receiver/GCS loss tests; preserve FS_GCS_ENABL=1')
add('Q_RTL_MODE',1,'Configure normal RTL to finish with VTOL return/landing instead of prolonged fixed-wing circling.','flying','Return path validated with real aircraft geometry and propulsion')
add('BATT_FS_LOW_ACT',6,'Replace warning-only with return, loiter down and QLAND.','bsource','Voltage calibrated and QRTL/return reserve demonstrated; these actions do not correct a bad current sensor')
add('BATT_FS_CRT_ACT',4,'Restore QLAND action that was disabled inside log17.','bsource','Voltage threshold calibration and sufficiently clear emergency landing area')
add('BATT_ARM_VOLT',46.8,'Restore the earlier minimum arming voltage from this same log; prevents treating the later 44 V allowance as a ready-to-fly standard.','bsource','12S chemistry and actual voltage calibration confirmed; this is not a state-of-charge estimate')

pending = {
 'ARSPD_TYPE':'Choose for the installed transducer; logged 0 means none. Do not copy SITL analog type 2.',
 'ARSPD_USE':'Propose 1 only after the physical sensor is healthy and calibrated. Absent from log17 because TYPE=0; SITL currently 0.',
 'ARSPD_RATIO':'Determine from physical pitot installation/calibration; do not copy SITL or peer value.',
 'SERVOx_FUNCTION':'Confirm direct vs inverse-V-tail mixing. Log17 has elevator19/rudder21 and all outputs9-16 disabled, so it does not establish a working CAN ruddervator setup.',
 'SERVOx_MIN/MAX/TRIM/REVERSED':'Measure actual travel and direction. Log1100-1900 != SITL1000-2000; simulated elevator1575 trim is not a hardware trim.',
 'TRIM_THROTTLE':'Logged45 retained as untuned placeholder; measure winged level flight. SITL requested25.4188 is stored as integer25 and uses a fictitious120 N pusher.',
 'PTCH_TRIM_DEG':'Retain logged0 until winged trim measurements; do not copy diagnostic -0.1485 deg.',
 'TECS_CLMB_MAX/TECS_SINK_MIN/TECS_SINK_MAX':'Measure sustainable climb/glide with real pusher and low battery; logged5/2/5 are unvalidated defaults.',
 'TECS_LAND_ARSPD':'Set slightly above the final validated AIRSPEED_MIN before an automatic landing approach; do not invent a speed while MIN is conditional.',
 'RTL_ALTITUDE/Q_RTL_ALT/WP_LOITER_RAD/RTL_RADIUS':'Location-dependent. Logged -1/15/60/0 does not define a reviewed first-forward-flight return path. At 25 m/s indicated and rho1.021, 25 deg bank gives ~164 m still-air radius. Start evaluating ~200 m loiter/return radius, 80 m return altitude and40 m VTOL return height in site-specific SITL, not as unconditional hardware values.',
 'INS_HNTCH_*':'Log16 spectrum ~75-90 Hz warrants design work; log17 notch still disabled. Do not copy GreatShark75 Hz or Mugin51 Hz, or ESC-RPM mode with non-telemetry KDE ESCs.',
 'Q_M_THST_HOVER':'Relearn with wings/pusher installed; retain last learned0.227437 as historical wingless value, not a prediction for winged25 kg.',
 'Q_M_BAT_CURR_MAX/BATT_WATT_MAX':'Remain0 until real propulsion-current measurement is calibrated; recent logged current is not total propulsion current.',
}

retained_reasons = {
 'Q_TRANSITION_MS':'Keep5000 ms; fade AFTER minimum airspeed, not total transition time.',
 'LEVEL_ROLL_LIMIT':'Keep5 degrees; level-transition bank cap after adding Q_OPTIONS bit0.',
 'Q_ASSIST_ANGLE':'Keep30 degrees; trigger becomes available with positive Q_ASSIST_SPEED.',
 'Q_ASSIST_DELAY':'Keep0.5 s; verify trigger timing in simulation.',
 'Q_BACKTRANS_MS':'Keep3000 ms; check return deceleration before changing.',
 'Q_TRANS_DECEL':'Keep2 m/s² as a tested-in-future stopping estimate; not an achieved measurement.',
 'Q_ENABLE':'Keep1 on hardware; simulation0 is a different fixed-wing-only setup.',
 'Q_FRAME_CLASS':'Keep1 QUAD; verify physical motor mapping.',
 'Q_FRAME_TYPE':'Keep1 X; verify motor order/direction.',
 'Q_M_THST_EXPO':'Keep0.75, doc starting value for20+inch props; no stand data justifies a new exponent.',
 'Q_M_PWM_MIN':'Keep1100; programmed KDE RANGE endpoint, still verify all ESCs.',
 'Q_M_PWM_MAX':'Keep1940; programmed KDE RANGE endpoint, still verify all ESCs.',
 'Q_M_PWM_TYPE':'Keep0 normal PWM; do not copy DShot settings.',
 'Q_RC_SPEED':'Keep400 Hz within KDE50-600 Hz input range.',
 'Q_M_BAT_VOLT_MAX':'Keep50.4 V for12S4.2 V/cell configuration.',
 'Q_M_BAT_VOLT_MIN':'Keep39.6 V compensation floor; not a battery-failsafe threshold.',
 'Q_M_BAT_IDX':'Keep0 first battery; GreatShark1 uses a different monitor topology.',
 'Q_M_HOVER_LEARN':'Keep2; relearn after weight changes.',
 'Q_M_SPIN_ARM':'Keep0.1 pending repeat bench confirmation; not inferred from motorKv.',
 'Q_M_SPIN_MIN':'Keep0.15 pending repeat bench confirmation.',
 'INS_GYRO_FILTER':'Keep20 Hz;20+inch initial recommendation.',
 'INS_ACCEL_FILTER':'Keep10 Hz; current documented VTOL starting range.',
 'ARMING_SKIPCHK':'Keep0; no bypass of checks.',
 'ARMING_REQUIRE':'Keep1.',
 'STALL_PREVENTION':'Keep1, but correct minimum airspeed still essential.',
 'SERVO_AUTO_TRIM':'Keep0 for the first diagnostic flight; controlled later trim procedure after measured surface setup.',
 'STAB_PITCH_DOWN':'Keep2 degrees on hardware; SITL0 was an intentional test simplification.',
 'FS_GCS_ENABL':'Keep1 and verify heartbeat failsafe; a disarmed QLAND event in log17 is not an airborne test.',
}
rows=[]
keys=set(proposal)|set(retained_reasons)|{'ARSPD_TYPE','ARSPD_USE','ARSPD_RATIO','Q_M_THST_HOVER','TRIM_THROTTLE','PTCH_TRIM_DEG','RTL_ALTITUDE','Q_RTL_ALT','WP_LOITER_RAD','RTL_RADIUS'}
for dictionary in [logged,sitl]:
    keys.update(k for k in dictionary if k.startswith(('Q_A_RAT_', 'Q_A_ACC_', 'Q_A_ANG_', 'RLL_RATE_', 'PTCH_RATE_', 'SERVO')))
keys.update(manifest['parameter_overlay'])
keys.update(k for k in logged if k.startswith(('TECS_','FLTMODE','BATT_FS','BATT_LOW','BATT_CRT')))
for k in sorted(keys):
    entry=proposal.get(k)
    status='PROPOSED_CHANGE' if entry and (k not in logged or abs(logged[k]-entry['value'])>1e-5) else ('EXPLICIT_RETAIN' if entry or k in retained_reasons else 'COMPARE_ONLY')
    if k in pending:status='MEASUREMENT_REQUIRED'
    reason=entry['reason'] if entry else retained_reasons.get(k,'')
    if k.startswith('Q_A_') and not reason:reason='Preserve established hover controller; first transition is not a reason to replace it with peer gains.'
    if k.startswith(('RLL_RATE_','PTCH_RATE_')) and not reason:reason='Initial hardware default only; verify reset defaults before SCALING_SPEED change; obtain a real fixed-wing tune. Do not substitute the diagnostic FDM gains.'
    if k.startswith('SERVO') and not reason:reason='Physical output/mixing/throw must be bench measured; current virtual FDM mapping is not a hardware calibration.'
    rows.append({'parameter':k,'log17_last_4_7_0':logged.get(k,'NOT_LOGGED'),'sitl_observed_4_7_1':sitl.get(k,'NOT_ACTIVE_OR_NOT_LOGGED'),'sitl_requested_overlay':manifest['parameter_overlay'].get(k,''),'candidate_4_7_1':entry['value'] if entry else 'RETAIN/PENDING','status':status,'why':reason or pending.get(k,'Retain for first FBWA; review before automatic flight.'),'condition':entry['condition'] if entry else pending.get(k,''),'greatshark_explicit':peers['Foxtech_GreatShark'].get(k,'NOT_SPECIFIED'),'mugin_explicit':peers['Mugin_EV350'].get(k,'NOT_SPECIFIED'),'source':entry['source'] if entry else ''})
with (root/'parameter-comparison.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
(root/'candidate.json').write_text(json.dumps({'status':'DRAFT_NOT_APPLIED_NOT_FLIGHT_VALIDATED','firmware':'Plane-4.7.1','assumed_mass_kg':25,'proposed':proposal,'pending':pending,'retained':retained_reasons},indent=2)+'\n')
(root/'sources.json').write_text(json.dumps({'review_date':'2026-10-09','official_sources':sources,'peer_parameter_files':identities},indent=2)+'\n')
(root/'fdm-evidence.json').write_text(json.dumps({'config':manifest['config'],'aircraft':manifest['aircraft'],'trim_control':manifest['trim_control'],'requested_parameters':manifest['parameter_overlay'],'firmware':manifest['firmware'],'limitations':manifest['limitations']},indent=2)+'\n')
(root/'candidate-review.param.txt').write_text('# READ-ONLY DISCUSSION DRAFT. All lines commented deliberately. Not an aircraft upload file.\n# Assumptions, pending measurements and return-path conditions: see README.md and candidate.json.\n'+''.join(f'# {k} {e["value"]}  # {e["reason"]}\n' for k,e in proposal.items()))
changed=[r for r in rows if r['status']=='PROPOSED_CHANGE']
header='| Parameter | Latest hardware log | Actual SITL | Proposed | Why / condition |\n|---|---:|---:|---:|---|\n'
table=header+''.join(f'| `{r["parameter"]}` | {r["log17_last_4_7_0"]:.6g} | {r["sitl_observed_4_7_1"]:.6g}' .replace('','')+' | '+str(r['candidate_4_7_1'])+' | '+r['why']+' '+r['condition']+' |\n' if isinstance(r['sitl_observed_4_7_1'],(int,float)) else f'| `{r["parameter"]}` | {r["log17_last_4_7_0"]:.6g} | inactive/unlogged | {r["candidate_4_7_1"]} | {r["why"]} {r["condition"]} |\n' for r in changed)
(root/'proposed-changes.md').write_text('# Proposed changes — none applied\n\n'+table)
calc={'mass_kg':25,'wing_area_m2':1.805,'span_m':3.95,'rho_kg_m3':1.021,'wing_loading_kg_m2':25/1.805,'aspect_ratio':3.95**2/1.805,'tas_to_eas':math.sqrt(1.021/1.225),'candidate_min_eas':18,'maximum_stall_eas_for_25pct_margin':18/1.25,'required_clmax_at_25kg':2*25*9.80665/(1.225*1.805*(18/1.25)**2),'cruise_tas_for_25_eas':25/math.sqrt(1.021/1.225),'turn_radius_at_25_eas_bank25':(25/math.sqrt(1.021/1.225))**2/(9.80665*math.tan(math.radians(25))),'stall_sensitivity':[{ 'assumed_clmax':cl,'stall_eas':math.sqrt(2*25*9.80665/(1.225*1.805*cl)), 'stall_tas':math.sqrt(2*25*9.80665/(1.021*1.805*cl))}for cl in [.8,1,1.2,1.5]]}
(root/'geometry-calculations.json').write_text(json.dumps(calc,indent=2)+'\n')
# Audit each proposed name against logged controller or observed 4.7.1 firmware;
# inactive Q groups checked explicitly in pinned source for enum/bitmask semantics.
missing=[k for k in proposal if k not in logged and k not in sitl]
assert not missing,missing
assert proposal['Q_OPTIONS']['value']==262177
assert proposal['AIRSPEED_MAX']['value']>=1.5*proposal['AIRSPEED_MIN']['value']
assert all(line.startswith('#') or not line.strip() for line in (root/'candidate-review.param.txt').read_text().splitlines())
print(json.dumps({'comparison_rows':len(rows),'proposed_numeric_changes':len(changed),'explicitly_retained_in_proposal':[k for k in proposal if k not in {r['parameter']for r in changed}],'geometry':calc},indent=2))
