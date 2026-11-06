# Project Spearhead — First Transition and Forward Flight

ArduPilot 4.7.1 parameter assessment and validation report

Reviewed 9 October 2026. Target firmware: **official ArduPlane 4.7.1**, source commit `dbe792162d06cab66c3475fd5556bf7a120f119e`.

**Subsequent control-data update:** a valid expanded export was received and integrated after this review snapshot. See [control-data-update.md](control-data-update.md) for the update and validation. The evidence and comparisons below retain their original snapshot; references to pending corrected controls should be read in that context.

**Recommendation:** preserve the established VTOL controller, add transition timeout and assistance, provide an explicit FBWA/VTOL recovery switch arrangement, and use a measured-airframe speed envelope. The numerical speed proposal below is conditional. Neither the hover logs nor the present FDM establishes the winged aircraft's stall speed or flight-ready CG. No aircraft, running simulator, GUI configuration or firmware parameter was changed by this review.

## Review summary

This report consolidates the proposed first-transition setup, every proposed parameter change and its rationale, comparisons against the newest supplied aircraft log and current FDM/SITL, two comparable aircraft, and the staged validation plan.

- **24 proposed numerical changes**, with conditions; no parameters have been applied.
- **291 relevant parameter comparisons**, backed by the complete extracted log and SITL parameter snapshots.
- Candidate airspeeds remain conditional on the winged aircraft’s measured weight, CG, airspeed calibration and a justified stall-speed bound.
- The present successful FDM tests cover the provisional reference aircraft in fixed-wing mode; a matching winged KDE QuadPlane transition remains to be validated.

Read the proposed-change table together with its conditions. The report describes a review candidate, not a completed aircraft tune. The folder’s README provides an index of every supporting deliverable.

## Evidence and configuration differences

| Item | Latest supplied aircraft log | Current successful FDM/SITL tests | First transition requirement |
|---|---|---|---|
| Evidence | `00000017.BIN`, supplied file dated 25 Sep; 1,241 distinct parameters | 9 Oct `roll_plus_28` manifest and its actual DataFlash parameters | Fresh board export after upgrading to 4.7.1 |
| Firmware | ArduPlane V4.7.0, `1511f271` | ArduPlane V4.7.1, `dbe79216` | Use the agreed official 4.7.1 on both |
| Aircraft | Quad-X hover configuration; no forward-flight mode recorded | `reference_design`, provisional 25 kg, reference MAD lift hardware, 120 N ideal pusher | Winged PT1 with KDE lift hardware and actual forward propulsion |
| QuadPlane | `Q_ENABLE=1` | `Q_ENABLE=0` | `Q_ENABLE=1`; transition still needs simulation validation |
| Airspeed sensor | `ARSPD_TYPE=0`, second sensor also 0; no ARSP messages | Analog simulated sensor type 2, **ARSPD_USE=0** | Choose actual sensor type; calibrate before enabling use |
| Aero model | Not identified by firmware log | Only `tail_update_2026_10_07`; rho 1.021; source moment centre (0,0,0) | Corrected extended control data and measured geometry |
| Mass / inertia | Log does not measure AUW / inertia | 25 kg; diag(2.5,4,5) kg m², provisional | Measure winged AUW, CG and inertia estimate |
| CG | No physical CG measurement in log | Diagnostic (0.013596,0,0.45) m | Do not use this as an aircraft CG recommendation |
| Surface mapping | Functions 4/19/70/21 on outputs1–4; outputs9–16 disabled | Virtual separate aileron/elevator/rudder; assumed ±20°, 1000–2000 µs | Verify inverse-V-tail/ruddervator mixing, CAN outputs, reversals and actual travel |
| Forward propulsion | A throttle output function does not prove a pusher was installed | Ideal linear 120 N maximum thrust | Installed engine/motor/propeller curve and response |

The `pt1_as_built.yaml` profile has the specified KDE5215XF-220, UAS75HVC+ and 21.5×7.3 propellers, but **pusher.installed=false**. The successful fixed-wing tests use the reference profile, not this PT1 configuration. The reference profile has no usable lift performance curve for running its MAD motors. Consequently, merely setting `Q_ENABLE=1` cannot turn the existing fixed-wing validation into a validated PT1 transition.

The simulation's source loads are applied at the user-specified moment centre (0,0,0). The base CSV header disagrees about its reference point; that remains a data clarification. Changing the moment-reference point is a coordinate transformation, not a way to change the aircraft's physical CG or make an unstable aircraft stable.

## What the newest hardware log actually supports

The log begins with the aircraft already armed around boot 207 s; the first disarm is 346.42 s. The later recorded arm/disarm is 1320.62–1387.22 s. Modes are QSTABILIZE/QHOVER/QLOITER, then a **disarmed** GCS failsafe switches to QLAND at 1580.60 s. No FBWA, CRUISE or transition was flown. A log snapshot cannot establish what is on the board today.

For a bounded late-hover window at boot 1335–1375 s, roll/pitch/yaw absolute tracking-error medians are **0.542° / 0.123° / 1.097°**, with 95th percentiles 1.668° /0.390° /3.364°. Motor-output medians are 1642/1665/1585/1581 µs. Yaw output still has positive bias, median 0.0867, but this is much less separated than the earlier log 16 example. This is enough to preserve hover gains as the starting point; it does not validate the heavier winged configuration. These statistics include pilot manoeuvres and are not a formal hands-off hover test.

At boot 1207.14 s, **BATT_FS_CRT_ACT changes 4→0**. Minimum arming voltage is also lowered from 46.8 to 44 V during this log. The last-hover window reaches42.49 V, with a median 42.75 V. Recorded current is only about 0.175 A median, inconsistent with total four-motor propulsion current. Capacity and current-based reserve estimates cannot be used from this monitor configuration without correcting the measurement. All first/last values and within-log changes are preserved in the accompanying JSON files.

## Geometry, weight and speed proposal

Repository wing references: **span 3.95 m, area 1.805 m², reference chord 0.457 m**. At the provisional 25 kg, wing loading is **13.85 kg/m²**, aspect ratio 8.64. The updated tail chord/area are still unset; the old 0.28 m tail chord is not substituted.

Use indicated/equivalent speed for parameter planning. For a calibrated pitot and rho 1.021, the approximate relation is `EAS = TAS × sqrt(1.021/1.225) = 0.912945 × TAS`. Thus 28 m/s FDM true speed corresponds to **25.56 m/s equivalent**, and a 25 m/s indicated cruise candidate corresponds to **27.38 m/s true**. ArduPilot's parameter description calls cruise indicated/calibrated/apparent airspeed. The simulator presently uses ArduPilot's standard atmosphere independently of the FDM's fixed density, so its sensor indication must be reconciled before transferring these numbers. [4.7.1 parameter source](https://github.com/ArduPilot/ardupilot/blob/dbe792162d06cab66c3475fd5556bf7a120f119e/ArduPlane/Parameters.cpp), [sensor calibration](https://ardupilot.org/plane/docs/calibrating-an-airspeed-sensor.html).

**Conditional candidate:** MIN 18, ASSIST 15, CRUISE 25, MAX 27 m/s indicated/equivalent. These values require the winged aircraft's stall speed to be **no higher than 14.4 m/s indicated**, so MIN has a 25% margin and ASSIST remains above stall. At 25 kg that requires CLmax at least **1.069**. We do not have a validated CLmax. The docs'25% margin is used here rather than the 20% phrasing in 4.7.1 parameter metadata. [Airspeed setup](https://ardupilot.org/plane/docs/airspeed-parameters-setup.html), [QuadPlane operation guidance](https://ardupilot.org/plane/docs/quadplane-reliability.html).

| Assumed CLmax — sensitivity only | Stall EAS | Stall TAS at rho 1.021 | Minimum at 1.25×stall EAS |
|---:|---:|---:|---:|
|0.8|16.65|18.24|20.81|
|1.0|14.89|16.31|18.61|
|1.2|13.59|14.89|16.99|
|1.5|12.16|13.32|15.20|

For example, CLmax 1.0 would already invalidate MIN 18. If weight changes, stall speed scales with sqrt(weight), provided configuration and CLmax stay the same. Recompute `MIN >=1.25×Vs`, and choose `ASSIST≈MIN−3` only if it remains above Vs; do not simply scale a peer's minimum speed by weight alone.

MAX 27 is an **automatic demanded-speed ceiling**, not a guarantee against overspeed in FBWA or a request to fly 27 indicated. It meets the 4.7.1 metadata guideline MAX≥1.5×MIN for this example, and still requires evidence of sustainable level-flight capability. For current validation, stay at or below **28 m/s FDM true speed**. No new simulation at 30 m/s or above was run. First steady flight would target roughly 25 indicated/27.4 true, conditional on measured airspeed calibration and stall margin.

## Every proposed numerical change

These are changes **relative to the final values in log 17**, not changes already applied. The complete row-level evidence and conditions are in [proposed-changes.md](proposed-changes.md), with actual SITL values alongside. [parameter-comparison.csv](parameter-comparison.csv) compares 291 relevant entries, including all current runtime-overlay entries, rate gains and servo settings. `NOT_ACTIVE_OR_NOT_LOGGED` must not be read as a numerical zero or an inferred default.

| Parameter | Log17 → proposed 4.7.1 | Reason |
|---|---|---|
|AIRSPEED_MIN|9→18|Conditional stall-margin example above; wingless9 is unsupported for winged flight.|
|Q_ASSIST_SPEED|−1→15|Enable low-speed and other assist triggers outside transitions; threshold above assumed stall.|
|AIRSPEED_CRUISE|12→25|Target near the tested28 true-speed point after density conversion.|
|AIRSPEED_MAX|22→27|Conditional auto-demand ceiling, consistent with MIN 18; not a flight-test target.|
|SCALING_SPEED|15→25|Choose near cruise before a new fixed-wing tune, with rate gains verified at factory defaults.|
|Q_TRANS_FAIL|0→20 s|Bound failure to reach transition airspeed; verify acceleration fits inside20 s.|
|Q_TRANS_FAIL_ACT|0→1|QRTL rather than QLAND on timeout, conditional on a rehearsed return path.|
|Q_OPTIONS|262144→262177|Preserve bit18 VTOL-only arming; add bit0 level transition and bit5 VTOL RC-failsafe QRTL.|
|ROLL_LIMIT_DEG|45→25|Gentle first-flight bank demands; still-air turn radius about 164 m at the candidate cruise.|
|PTCH_LIM_MAX_DEG|20→10|Constrain initial nose-up demand pending climb capability.|
|PTCH_LIM_MIN_DEG|−25→−10|Constrain initial nose-down demand.|
|PTCH2SRV_TCONST|0.5→0.75 s|Slower initial pitch-demand shaping; both peer reference files use0.75.|
|RLL2SRV_RMAX|0→45 deg/s|Bound initial roll-rate demand.|
|PTCH2SRV_RMAX_UP|0→22 deg/s|Initial pitch-rate limit based on560/25 guidance.|
|PTCH2SRV_RMAX_DN|0→22 deg/s|Same downward demand limit.|
|FLTMODE3|18→19|QLOITER recovery slot before FBWA.|
|FLTMODE4|18→5|FBWA for first manual-forward-throttle transition.|
|FLTMODE5|19→21|Explicit QRTL return slot.|
|FLTMODE6|19→20|Explicit QLAND slot.|
|FS_LONG_ACTN|0→1|Explicit RTL on long failsafe even from automatic modes.|
|Q_RTL_MODE|0→1|Make RTL complete with VTOL return/landing.|
|BATT_FS_LOW_ACT|0→6|Return/loiter down/QLAND rather than warning-only, after voltage/reserve validation.|
|BATT_FS_CRT_ACT|0→4|Restore critical-battery QLAND action disabled in log 17.|
|BATT_ARM_VOLT|44→46.8 V|Restore earlier12S minimum arming voltage from the same log, after calibration/chemistry confirmation.|

The transition timer and options were cross-checked in the **pinned4.7.1 source**, not only today's wiki. `Q_TRANSITION_MS=5000` is retained: it is the lift fade after reaching minimum airspeed. The timeout applies to waiting for airspeed. Bit0 restricts bank to `LEVEL_ROLL_LIMIT=5` and changes how climb is handled during transition. Bit19 is not added. [Transition behaviour](https://ardupilot.org/plane/docs/quadplane-transitions.html), [4.7.1 options and action enums](https://github.com/ArduPilot/ardupilot/blob/dbe792162d06cab66c3475fd5556bf7a120f119e/ArduPlane/quadplane.cpp).

QRTL is **not guaranteed to be a pure VTOL abort**: its default distant return can use fixed-wing flight. An aborted transition with inadequate forward thrust must not be followed by an untested automatic transition attempt. Validate near/far return behaviour, or deliberately choose a rehearsed QLAND/pure-VTOL return strategy for the test site's available area and energy. The QRTL choices above remain conditional until that review. [QRTL behaviour](https://ardupilot.org/plane/docs/qrtl-mode.html).

Mode slots1/2 stay QSTABILIZE17/QHOVER18, RC channel7 stays the mode switch, and arming checks stay enabled. Verify these six positions on the actual radio. The battery action numbers are Plane-specific; they are not Copter action enums. [4.7.1 battery action definitions](https://github.com/ArduPilot/ardupilot/blob/dbe792162d06cab66c3475fd5556bf7a120f119e/libraries/AP_BattMonitor/AP_BattMonitor_Params.cpp).

## What I would preserve or measure rather than copy

Keep hover roll and pitch P/I/D **0.25/0.25/0.0036**, yaw P/I/D **0.5/0.018/0**, gyro/accel filters 20/10 Hz, frame Quad-X, normal PWM, Q_RC_SPEED400 Hz, KDE endpoints 1100/1940, Q_M_THST_EXPO0.75, voltage compensation50.4/39.6 V and battery index 0. Reconfirm motor-start endpoints and relearn hover thrust with the final weight. The21.5-inch props fit the documented20+inch initial exponent/filter guidance. Q_M_THST_HOVER0.227437 is a normalized thrust-control value, **not** the manufacturer test sheet's throttle percentage. [VTOL tuning process](https://ardupilot.org/plane/docs/quadplane-vtol-tuning-process.html), [KDE ESC input range](https://www.kdedirect.com/products/kde-uas75hvc).

The hardware fixed-wing gains are essentially factory starting values, without logged fixed-wing validation. Keep/reset those starting gains for the first new tune; do not load the FDM gains. Choosing SCALING_SPEED25 changes the effective tune and therefore requires default-rate-gain verification and a fresh physical tune. Set it before tuning. [Speed scaling](https://ardupilot.org/plane/docs/airspeed-parameters-setup.html).

| Fixed-wing parameter | Hardware log 17 | Actual SITL 4.7.1 | Proposed physical-aircraft starting point |
|---|---:|---:|---|
|RLL_RATE_FF|0.345|0.08625|Verified4.7.1 default 0.345, then tune|
|RLL_RATE_P|0.08|0.02|Default 0.08, then tune|
|RLL_RATE_I|0.15|0.075|Default 0.15, then tune|
|RLL_RATE_D|0|0|Default 0, then tune|
|PTCH_RATE_FF|0.345|0.45|Default 0.345, then tune|
|PTCH_RATE_P|0.04|0.08|Default 0.04, then tune|
|PTCH_RATE_I|0.15|0.15|Default 0.15, then tune|
|PTCH_RATE_D|0|0|Default 0, then tune|
|PTCH2SRV_TCONST|0.5|0.2|Candidate0.75; validate response|
|SERVO2_TRIM|1500|1575|Measure actual ruddervator trim|
|TRIM_THROTTLE|45|**25**, requested25.4188|Measure level-flight throttle|
|PTCH_TRIM_DEG|0|−0.148522|Measure actual winged level trim|

The SITL trim-throttle parameter is integer-valued in this firmware: its logged stored value is 25, despite the manifest request25.4188. The full comparison preserves both. Its gain scaling depends on assumed±20° virtual travel. Physical1100–1900 endpoints, unknown mechanical travel, inverse-V-tail mixing and a different CG can substantially change control authority. A simulation recovery pass cannot choose hardware gains without those measurements. [Current fixed-wing tuning procedure](https://ardupilot.org/plane/docs/new-roll-and-pitch-tuning.html).

Leave exact sensor model/ratio, ruddervator functions/reversals/travel, real cruise throttle, climb/sink capabilities, automatic approach speeds and site return altitudes pending measurement. The current60 m loiter radius is too small for an unqualified25° bank/25 indicated return assumption: nominal still-air radius is 164 m, before wind and tracking margin. Evaluate roughly 200 m radius and suitable obstacle-clearing altitudes in site-specific simulation; no location-independent altitude or radius has been applied.

Do not invent a notch from prop diameter. Earlier log 16 batch data had75–90 Hz peaks; log 17 still has the notch disabled. Use the actual final configuration's spectra across throttle levels. KDE HVC ESCs provide no ESC-RPM telemetry in our setup. Current, capacity and current limiting need a calibrated propulsion-current measurement first.

## Similar aircraft with published parameter sets

These are upstream reference files, not complete exports from the exact flights shown by the manufacturers. Their metadata does not establish flight weight, firmware version or acceptance testing for each tune. Copies and hashes from the locally pinned4.7.1 checkout are saved here; comparisons use only explicitly listed values, not inferred peer defaults.

| Aircraft | Relevant similarity | Important difference |
|---|---|---|
|Foxtech GreatShark330|Manufacturer lists3.2 m span, about 24 kg MTOW, four 22-inch lift props, tail electric pusher and airspeed sensor; upstream file uses V-tail functions79/80.|Our wing is 3.95 m and provisional 25 kg, KDE21.5-inch propulsion; same wing loading and tail effectiveness are not established.|
|Mugin EV350|Manufacturer currently lists3.5 m span,23 kg MTOW, four electric lift motors and rear electric cruise propulsion.|Wing area 1.0 m² vs our1.805; published stall16 m/s and cruise21–24 m/s have unspecified reference conditions. Current2026 product may differ from the older upstream tune.|

[GreatShark manufacturer](https://www.foxtechfpv.com/foxtech-great-shark-330-vtol.html), [Mugin manufacturer](https://www.muginuav.com/product/mugin-ev350-carbon-fiber-full-electric-vtol-uav-platform/).

| Explicit parameter | GreatShark reference | EV350 reference | Our candidate / disposition |
|---|---:|---:|---|
|AIRSPEED_MIN|18|18|18 conditional on our own stall margin|
|Q_ASSIST_SPEED|15|15|15 conditional|
|AIRSPEED_CRUISE|23|23|25 conditional near present FDM point|
|AIRSPEED_MAX|35|28|27 conditional demand cap|
|TRIM_THROTTLE|60|50|Measure ours; don't copy either|
|PTCH2SRV_TCONST|0.75|0.75|Candidate0.75|
|Q_TRANS_FAIL|20|not specified|20 conditional acceleration/abort test|
|Q_M_THST_HOVER|0.338|0.264|Wingless learned 0.227437; relearn winged|
|Q_M_BAT_IDX|1|0|Retain0 for our monitor topology|
|Q_M_THST_EXPO|not specified|0.7|Retain0.75 pending thrust stand evidence|
|INS_HNTCH_FREQ|75|51|Derive from our measured final configuration|
|PTCH_RATE_FF|1.835124|1.093860|Neither is evidence for our elevator authority|
|RLL_RATE_FF|0.442834|0.541178|Neither is evidence for our aileron authority|
|ROLL_LIMIT_DEG|40|40|25 for initial gentle manoeuvres|
|Q_RTL_ALT|40|40|Site-specific; logged 15 is not reviewed for our first transition|
|RTL_RADIUS|200|180|Evaluate~200 with the actual field and wind|

[GreatShark upstream parameters](https://github.com/ArduPilot/ardupilot/blob/master/Tools/Frame_params/QuadPlanes/Foxtech_GreatShark.param), [EV350 upstream parameters](https://github.com/ArduPilot/ardupilot/blob/master/Tools/Frame_params/QuadPlanes/Mugin_EV350.param). Their substantial gain differences despite similar overall size illustrate why geometry/weight alone cannot determine our gains.

## KDE propulsion implications

At 25 kg, each rotor must provide about 61.3 N in level hover before control reserve. The integrated manufacturer static table gives, for21.5×7.3 two-blade props, about 79.1 N maximum per motor at 42 V and95.88 N at 50.4 V. Four-rotor static thrust/weight is therefore about 1.29 and1.56 respectively, before installation, inflow, battery sag and environmental effects. The motor's180-second limits are 44 A and1955 W; the 50.4 V full-throttle row exceeds both at 55.8 A/2811 W. An ESC's75 A rating does not increase the motor's limit. This supports maintaining battery reserve and bounding unsuccessful transitions, not inventing higher motor limits. [KDE performance sheet](https://cdn.shopify.com/s/files/1/0496/8205/files/KDE5215-220kV.pdf?v=1680823195), [motor limits](https://www.kdedirect.com/products/kde5215xf-220).

The model has static table interpolation and an idealized forward motor, not a validated inflow/propwash model. Do not infer precise transition time, hover battery reserve or torque margin from these static thrust ratios. Weight, propeller clearance and final structure differ from the wingless logged flights.

## Validation sequence before the first transition

1. **Complete the physical configuration record.** Weigh the winged aircraft, measure CG in the AeroDB axes, identify forward propulsion, measured tail geometry and surface travel. Resolve the pending corrected control data. Establish a justified stall-speed bound for the first test configuration; don't start by probing stall in an unvalidated transition.
2. **Verify4.7.1 and bench mapping.** Obtain a fresh full parameter export. Check that calibrated sensors/RC/board orientation survive the upgrade. Verify all motor directions/order, pusher idle/full command, stabilizing surface direction, inverse-V-tail mixing and mode slots with motors made safe. Set physical pitot type; verify health, zero and pressure response. Keep sensor use disabled until calibration checks are satisfied.
3. **Build the matching QuadPlane FDM configuration.** Use KDE lift curves, measured AUW/CG, actual forward propulsion, calibrated surface mapping and consistent indicated/true speed. Keep the existing fixed-wing test evidence unchanged. Hover and transition tests must use Q_ENABLE1. Initial validation stays within the requested≤28 m/s true-speed scope.
4. **Exercise transitions and recovery in simulation.** Test QHOVER/QLOITER→FBWA, reaching the validated MIN before timeout, the subsequent5-second fade and FBWA→QHOVER/QLOITER. Record altitude loss, airspeed, rate/attitude tracking, thrust/PWM saturation and control-data extrapolation. Inject inadequate pusher thrust, misleading/unhealthy airspeed, representative wind, battery sag, RC/GCS loss and near/far QRTL. Verify the intended abort mode without an unintended new transition. A stock QuadPlane SITL test can check firmware logic but does not validate Spearhead aerodynamic performance.
5. **Recheck hover with wings installed.** Confirm controller tracking, motor reserve, yaw bias and vibration at the actual mass. Relearn hover thrust. Do not substitute the lighter wingless learned value as a 25 kg prediction.
6. **First forward segment: pilot-commanded FBWA.** Use a field and height that allow the rehearsed VTOL recovery. Maintain a into-wind, nearly level transition with measured airspeed called out; pilot controls forward throttle in FBWA. Keep manoeuvres gentle, observe the validated speed band and verify lift fade. The configured cruise/MAX values do not enforce speed in this mode. Rehearse stick behaviour on returning to QHOVER, where throttle becomes vertical-speed demand. Exact altitude, corridor, battery reserve and stop criteria belong to the tested site plan. [Flying a QuadPlane](https://ardupilot.org/plane/docs/quadplane-flying.html).
7. **Review the first log before tuning further.** Inspect CTUN.As and physical ARSP, transition messages, ATT/RATE/PIDP/PIDR, AETR surface demand and scaling, lift/pusher outputs, battery voltage/current and EKF/vibration. Confirm trim and actual rate response. Tune the pitch/roll loops using the current4.2+ procedure, then identify sustainable climb, sink and cruise throttle before adding CRUISE/TECS/AUTO and automatic approaches. [Current tuning](https://ardupilot.org/plane/docs/new-roll-and-pitch-tuning.html), [TECS guide](https://ardupilot.org/plane/docs/tecs-total-energy-control-system-for-speed-height-tuning-guide.html).

Current validation status: **research and parameter-name/value comparison complete; aircraft transition and proposed tune not tested.** Existing fixed-wing 28 m/s disturbance-recovery tests remain valid only for their recorded diagnostic configuration. Lower-speed failures to find trim in the partial control-data coverage are not measured stalls.

## Deliverables and audit

- `parameter-comparison.csv`:291 relevant parameters, latest logged value, actual SITL value, requested SITL overlay, candidate, peer explicit values and rationale.
- `proposed-changes.md`:24 numerical proposals with row-level conditions.
- `candidate.json`:proposals, retained settings and unresolved measurement-dependent fields.
- `candidate-review.param.txt`:all parameter lines commented; deliberately not an aircraft upload file.
- `log17-first-parameters.json`, `log17-last-parameters.json`, `log17-parameter-changes.json`, `log17-last-baseline.param`:historical log evidence; baseline is not a deployment overlay.
- `log17-summary.json`, `log17-last-hover-window.json`:message/mode/arm and bounded-hover statistics.
- `sitl-last-parameters.json`, `sitl-log-identity.json`, `fdm-evidence.json`:actual 4.7.1 firmware data and matching FDM manifest excerpt.
- `geometry-calculations.json`, `sources.json`, upstream peer parameter copies:reproducible calculations and provenance.

Validation performed: both DataFlash files parsed; runtime overlay compared with actual stored values; all 24 proposed names exist in logged/4.7.1 parameter evidence; Q_OPTIONS bits and action enums checked against pinned4.7.1 source; candidate MIN/MAX relationship and commented-draft format checked. No flight simulation was launched for this review.
