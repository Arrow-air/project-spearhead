---
title: "Spearhead PT1 and PT1.5 Flight Log Review, October 18, 2026"
sidebar_label: "Flight Log Review 2026-10-18"
sidebar_position: 9
---

# Spearhead PT1 and PT1.5 Flight Log Review

**Date:** October 18, 2026
**Author:** errrks (log analysis with Claude)
**Logs:** 00000014 (September 6), 00000015, 00000016, 00000017 (September 19), 00000033 (October 18). Pixhawk 6C, ArduPlane V4.7.0 through log 17, V4.7.2 in log 33. Same field for all five: 39.8582 N, 32.5461 E, 1,123 m elevation, south of Ankara. Local time is UTC plus 3. Times in the text are log seconds unless marked.
**Status:** Review. Findings are hypotheses with a stated confidence and a way to confirm. Nothing here has been verified on the aircraft yet.

Companion to SPH-E-006 (PT1 as-built note, open items AB-03, AB-08 and AB-12) and SPH-E-003 (parts list).

## 1. What the logs say about October 18

The hover lasted 19 s from arm to disarm and 9 s airborne. It ended with the aircraft sinking from 7.5 m at up to 2.3 m/s while the pilot held full climb stick, and a left yaw of about 100° that the controller could not stop. The sequence in the log:

| Log s | Event | Evidence |
|---|---|---|
| 554.9 | Armed in QLOITER, nose heading 84° (east). Pack 48.84 V at rest | BAT, ATT, EV 10 |
| 555.5 | Motors at idle, 1184 µs. Pack drops to 47.95 V | 0.86 V drop at idle. September idle drop at the same command was 0.22 to 0.26 V |
| 560.5 to 563.0 | Throttle stick 1168 to 1940 µs, motors spool, pack 47.0 to 43.5 V | RCIN C3, RCOU C5 to C8 |
| 563.5 | Liftoff. Pack 40.54 V. Motor 2 (rear left, CCW) at the 1898 µs ceiling from this moment on. Heading 83° | RFND 0.13 to 0.30 m, RCOU C6 1898 |
| 564.5 | "Weathervane Active: nose in". Heading already 78° and drifting left | MSG, ATT |
| 565.0 | Pack minimum 39.83 V | BAT |
| 565 to 569.5 | Climb to 7.54 m at up to 1.7 m/s. CCW minus CW motor split grows from 130 to 316 µs. Heading 73° to 46° | QTUN, RCOU |
| 570.3 | "Potential VTOL Thrust Loss (1)". Motor 1 at 1895, motor 2 at 1898, motor 3 at 1439, motor 4 at 1625 µs | MSG, ERR 25/1, RCOU |
| 570.5 to 573.0 | Throttle output at its limit (QTUN ThO 1.0). Split 420 to 437 µs. Descent 0.9 to 2.3 m/s. Stick 1870 to 2011 µs. Heading 35° to 355° | QTUN, RCOU, RCIN, ATT |
| 573.3 | "Battery 1 is low 41.04V used 42 mAh" (10 s under the 42 V threshold) | MSG |
| 573.5 | Touchdown at roughly 1 to 2 m/s. Arm switch (RC5) to low. Heading 347° | RFND 0.11 m, RCIN C5 |
| 573.8 | Disarmed. Pack 46.6 V at 0.2 s, 47.6 V at 8 s, 48.1 V at 36 s | EV 11, BAT |

The pusher motor was never commanded (output 3 stayed at 1100 µs for the whole log). The 70 m/s² lateral vibration peak in this log is on the ground at 561 to 562 s during spool with the gear loaded, not in flight. Airborne vibration was 6 to 8 m/s² on all axes.

## 2. Headline findings, ranked

| # | Finding | Confidence | How to confirm or act |
|---|---|---|---|
| 1 | **The hover ended because yaw and collective thrust were competing for the same two motors and the CCW pair ran out of range.** From liftoff, motor 2 sat at the 1898 µs ceiling and motor 1 joined it at 568 s. The yaw controller was asking for right yaw at 0.4 to 0.65 of full output the whole time and the nose still turned left at 6 to 25°/s. ArduPilot's mixer protects roll, pitch and yaw before collective, so it pulled motors 3 and 4 down to 1300 to 1500 µs to keep the differential, the average thrust fell below hover, and the aircraft sank with the stick at full. The thrust loss message is this saturation being detected on motor 1, not an ESC fault | High (RCOU, RATE.YOut, QTUN.ThO, ATT all agree) | None needed on the mechanism. The fixes are items 2, 3 and 7 |
| 2 | **The yaw torque bias is not new. It has been there since the first hover on September 6.** In every airborne segment the CCW pair (motors 1 and 2) runs 200 to 283 µs above the CW pair (motors 3 and 4), and the yaw output averages +0.21 to +0.31 of full authority in calm air with no pilot yaw input. The heading sits left of target on average in every flight (mean yaw error −3° to −18°). It was present with no wings (logs 14 and 16), so it is not aerodynamic. It is not a reversed motor either: one reversed motor would need the full differential and the aircraft could never hold a heading, but it did. The Quiver precedent is a 1.5 to 1.7° arm twist producing a +135 µs split, cured to +19 µs by shimming. Spearhead's split is twice that, so the equivalent misalignment is 2 to 3°, static or under load. The split varied between flights (+76 µs in the second flight of log 17, +283 in log 14), which points at something that moves: boom or spar flex under thrust, a mount that shifts, or prop tracking, rather than a fixed machining error | High on the observation, medium on motor tilt as the cause | Digital inclinometer on each motor mount about the boom axis, static first, then with the aircraft strapped down at hover thrust. Confirm each prop's handedness and rotation per position (1 front right CCW, 2 rear left CCW, 3 front left CW, 4 rear right CW). Shim until the airborne CCW minus CW split is under 50 µs. This is AB-12 |
| 3 | **The pack delivered about 3.5 times the voltage drop per amp compared with September.** The cleanest like for like number is the idle drop at the same 1184 µs command on the same props: 0.86 V on October 18 against 0.22 to 0.26 V on September 19. The hover sag was 8.3 V from rest (48.8 to 40.5 V) against 3 V in September (47.4 to 44.5 V), with higher current today. Recovery after disarm was slow (46.6 V at once, 48.1 V after 36 s), which is cold cell polarisation behaviour. The field was 8.6 to 9.1 °C and the airspeed probe read 9.7 °C. A cold pack is the simplest explanation, but a degraded cell group or a high resistance joint (QS8 socket, 8 AWG solder cups, FCHUB pads) cannot be separated from it without a current measurement | High that the source resistance was far higher than September, medium on cold as the only cause | Measure pack temperature before arming. Warm the pack to 25 to 30 °C before flight (insulated box, charge immediately before, or warm it in the car). Cell checker before and after each hover (SPH-E-003 line 20). Bench: pack resistance at 10 °C and 25 °C with a known load, and millivolt drop across each HV joint at 50 A |
| 4 | **The flight controller has no current measurement, and has not had one since September 6.** BATT1 reads 0.1 to 0.6 A throughout every log, including hover. BATT_MONITOR is 4 (analog) on the PM02 V3, and per the SPH-E-003 battery interface build the PM02 carries only the FC load from a pigtail, so its shunt never sees motor current. The FCHUB 440 A sensor on POWER2 is not configured (BATT2_MONITOR 0, AB-08 still open). Consequences: the mAh counter is meaningless (42 mAh for a hover), the resting voltage estimator reports 4 to 9 Ω, capacity failsafes are dead, and nobody can tell sag from a weak pack | High | Close AB-08. FCHUB V/I into POWER2, make it the primary monitor, calibrate BATT_AMP_PERVLT against a clamp meter, keep the PM02 as the FC supply only |
| 5 | **The battery failsafe thresholds have been lowered three times to get past warnings, and today's critical threshold was 3.5 s from commanding QLAND in the hover.** Arm voltage went 46.8 to 44.0 (changed in the field at log 17 1303 s) to 44.4 V. Low went 44.4 to 42.0 V. Critical went 42.0 to 40.8 V. The source is raw voltage (BATT_FS_VOLTSRC 0) with a 10 s timer. Today the pack was continuously below 40.8 V from 563.5 to 570.0 s, 6.5 s of the 10 s timer, then again from 572.5 s. With this pack at this temperature any hover longer than 10 s would have handed the aircraft to QLAND at what was effectively 80 % state of charge | High | Stop lowering the thresholds. Once current works, use the sag compensated source (BATT_FS_VOLTSRC 1) or set thresholds from a resting voltage policy plus a ground cell check. Decide the critical action deliberately, QLAND is reasonable |
| 6 | **The pack starts every session at 4.05 to 4.07 V per cell, not 4.20, and September 19 looked like a pack with far less than 22 Ah in it.** Starts: 48.61 V (September 6), 48.51 V (September 19 first log), 48.84 V (October 18). On September 19 about 5.5 minutes of total airborne time took the pack from 48.5 V to under 46.8 V resting (log 16 end) and then to about 44 V resting (3.67 V per cell) after log 17 flight 1, where it sat for 15 minutes failing the arming check. A healthy 22 Ah pack at the hover currents implied by the KDE table would have used about a third of its capacity and rested near 45.5 V. Either the charger terminates low, or the usable capacity is well below 22 Ah, or the hover current is much higher than the motor table suggests | Medium (no current data, voltage only) | Check the charger termination and balance state. Bench discharge the pack at 20 to 30 A to the vendor cutoff and record delivered Ah. Ask motorobit or ProFuse for the cell's rated max voltage |
| 7 | **Wind and the weathervane added yaw demand on top of the bias.** Reanalysis wind at the field at 17:00 local was 4.2 to 4.7 m/s from 249° with gusts to 9.7 m/s. The nose was at 84° at liftoff, so the aircraft hovered tail to wind with a 3.95 m wing and a V tail, which is the unstable direction for a weathercock moment. September 19 had 1.5 to 2 m/s. Q_WVANE_GAIN was changed from 0 to 1 for this flight and the weathervane engaged at 564.5 s, one second after the left drift had already started, so it did not start the departure but it was commanding yaw during it | Medium on the magnitude, high that the conditions differed from September | Hover tests nose into wind. Keep Q_WVANE_GAIN at 0 until the yaw bias is fixed and a calm hover shows a split under 50 µs. Pilot brief: wind limit for hover tests until then |
| 8 | **Thrust margin at 40 V is thin even with the yaw bias fixed.** Q_M_SPIN_MAX 0.95 puts the ceiling at 1898 µs. From the KDE 5215XF-220 sheet with the 21.5 x 7.3 DP prop, the four motors can give about 37.5 kgf at 50.4 V, 33.4 kgf at 44.5 V (the September hover bus) and 29.8 kgf at 40.5 V (today's). Today's airborne average was 1738 µs against 1546 to 1615 µs in September. The voltage drop alone costs about 55 µs per motor at constant thrust, the pusher pack and motor added roughly 2.7 kg, and the rest is the yaw fight and the wind. The aircraft has not learned its hover throttle since September (Q_M_THST_HOVER 0.227 unchanged) so the altitude loop feedforward was also low | Medium (KDE table, as flown mass not recorded) | Record the as flown mass for each configuration. Target a hover at or below 1650 µs average on a warm pack. If the mass stays at the PT1.5 figure, the 3 blade prop or a 14S bus are the levers, and both are PT2 questions |
| 9 | **Landing was hard and the airframe needs a look.** Descent 2.0 m/s at 573.0 s, 0.3 m above ground at 573.5 s, roll swing to −8.8° on touchdown, disarm by switch at touchdown (good, no prop strike while spinning) | High | Inspect gear, boom roots, motor mounts and prop tips before the next arm. A bent mount would also change item 2 |
| 10 | **Field parameter changes worth knowing.** Today: RC8_OPTION 31 to 0 at 87 s (motor emergency stop removed after a false trip at 65 s, this was the only HV kill per decision D2), FLTMODE1 and 2 to QLOITER, FLTMODE3 to FBWA, FLTMODE6 to QRTL. Since log 17: FS_GCS_ENABL 0 (after a GCS failsafe put the aircraft into QLAND on the ground at log 17 1580 s), FS_LONG_ACTN 1, INS_HNTCH enabled at 80 Hz throttle based, Q_ASSIST_SPEED and ALT 15, Q_RTL_MODE 1, Q_TRANS_FAIL 20 s, AIRSPEED_MIN 18 / CRUISE 20 / MAX 24, firmware 4.7.2. September 19: RC2_REVERSED set to 1 between flights in log 16, RTL_ALTITUDE 100 to 5 to −1, FENCE_ENABLE toggled | High | Re-enable the E-stop switch on RC8 before the next flight, or document the replacement kill path. Export the parameter file from the aircraft now so AB-03 can close against it |

## 3. Log inventory

| Log | Date, local | Firmware | Configuration | Armed segments | Airborne | Notes |
|---|---|---|---|---|---|---|
| 14 | Sep 6, 19:30 | 4.7.0 | No wings. KDE PWM set. Here4 as GPS1 (DroneCAN). RC protocol reported None at boot | 1 (log opens armed, LOG_DISARMED 0) | 46 s, 6.1 m | QSTABILIZE, QHOVER, QLOITER. First hover. Yaw error 22.8° rms, CCW minus CW +283 µs. Pack 48.6 to 44.7 V |
| 15 | Sep 19, 12:16 | 4.7.0 | Bench, USB only (BAT 1.3 V) | 0 | 0 | 14 s power up |
| 16 | Sep 19, 14:01 to 14:38 | 4.7.0 | No wings | 12 (ten short ground arms 94 to 732 s, two flights) | 62 s to 7.9 m (1412 to 1485 s), 95 s to 11.4 m (2032 to 2143 s) | Hover throttle learned 0.25 to 0.228. "Battery 1 is low 44.27V" at 2097 s. Arming check failed on voltage after landing |
| 17 | Sep 19, 15:41 to 16:08 | 4.7.0 | Wings and tail on | 2 (log opens armed) | 119 s to 17.7 m (224 to 343 s), 52 s to 12.7 m (1334 to 1386 s) | "Battery 1 is low 44.06V" at 244 s. Pack rested under 44 V for 15 min. BATT_FS_CRT_ACT set 0 at 1207 s and BATT_ARM_VOLT set 44 at 1303 s to arm again. GCS failsafe to QLAND on the ground at 1580 s |
| 33 | Oct 18, 17:18 to 17:29 | 4.7.2 | Wings, tail, pusher motor and 12S 15 Ah pusher pack fitted. F9P GPS1 serial, Here4 GPS2. MS4525 airspeed. LOG_DISARMED 1 | 3 (2.6 s at 206 s, 18.9 s at 555 s, 2.7 s at 610 s) | 9 s to 5.9 m | The hover in section 1. Motor E-stop tripped at 65 s and the option was removed at 87 s |

## 4. Per flight trend table

Airborne portions only (rangefinder above 0.5 m). Motors are ArduPilot numbers on outputs 5 to 8: 1 front right CCW, 2 rear left CCW, 3 front left CW, 4 rear right CW. YOut is the yaw mixer output, positive is nose right. Yaw error is actual minus target, negative is left of target.

| Log, flight | Config | Air s | Max AGL m | Pack V arm / min / air mean | M1 M2 M3 M4 mean µs | CCW − CW µs | M1 / M2 time at ceiling | YOut mean / max | Yaw err mean / rms ° | Yaw rate min / max °/s | ThO mean / max | Roll / pitch err rms ° | VIBE max x/y/z |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 14 | no wings | 46 | 6.1 | 48.61 / 44.71 / 45.23 | 1642 1655 1342 1389 | +283 | 0 / 0 | +0.31 / +0.70 | −17.7 / 22.8 | −27 / +23 | 0.21 / 0.28 | 1.5 / 1.8 | 18 / 8 / 16 |
| 16, f1 | no wings | 62 | 7.9 | 48.38 / 43.63 / 44.56 | 1669 1656 1389 1401 | +268 | 0 / 0 | +0.29 / +0.70 | −11.4 / 18.1 | −22 / +17 | 0.23 / 0.40 | 0.7 / 2.2 | 9 / 20 / 27 |
| 16, f2 | no wings | 95 | 11.4 | 48.01 / 43.84 / 44.47 | 1639 1639 1432 1445 | +201 | 0 / 0 | +0.21 / +0.46 | −2.6 / 7.4 | −25 / +14 | 0.23 / 0.34 | 0.4 / 0.7 | 8 / 17 / 25 |
| 17, f1 | wings | 119 | 17.7 | 47.40 / 43.01 / 43.51 | 1656 1637 1448 1441 | +201 | 0 / 0 | +0.23 / +1.16 (3.3 % above 0.9) | −5.1 / 15.1 | −10 / +28 | 0.23 / 0.37 | 0.7 / 1.2 | 14 / 16 / 26 |
| 17, f2 | wings | 52 | 12.7 | 46.14 / 42.35 / 42.78 | 1645 1662 1576 1577 | +76 | 0 / 0 | +0.09 / +0.18 | −1.2 / 1.6 | −3 / +7 | 0.28 / 0.40 | 0.8 / 0.4 | 7 / 5 / 23 |
| 33 | wings, pusher kit | 9 | 5.9 | 48.75 / 39.83 / 40.55 | 1862 1883 1537 1668 | +270, peak +437 | 43 % / 71 % | +0.40 / +0.65 | −27.0 / 30.5 | −25 / +1 | 0.77 / 1.00 | 1.7 / 2.5 | 15 / 70 (ground) / 18 |

Other channels in every flight: GPS 30 to 32 satellites, HDOP 0.47 to 0.56 (log 17 flight 1 dipped to 11 satellites on the Here4 alone). EKF variances well under the 0.8 threshold. CRSF link quality 99 to 100 %. Vcc 5.07 to 5.16 V. No accelerometer clipping.

## 5. Weather at the field

Open-Meteo at 39.8582 N 32.5461 E, model elevation 1,128 m. September from the ERA5 archive, October 18 from the forecast model's past day record. 10 m wind on a grid of several kilometres, so local gusts and direction can differ. The log barometer read 892.6 hPa today, within 2 hPa of the model, and the airspeed probe read 9.6 to 9.8 °C, so the model is representative.

| Session, local | Air temp °C | RH % | Surface pressure hPa | Wind m/s, from | Gust m/s | Cloud % | Air density kg/m³ |
|---|---|---|---|---|---|---|---|
| Sep 6, 19:30 | 24 to 26 | 25 to 30 | 894 | 5.3 to 6.1 from 354 to 003 (north) | 11.5 to 12.0 | 0 | 1.046 |
| Sep 19, 14:00 to 16:10 | 22.3 to 23.2 | 43 to 47 | 895.5 | 1.5 to 2.1 from 304 to 313 (northwest) | 5.4 to 5.7 | 27 to 100 | 1.055 |
| Oct 18, 17:18 | 8.6 to 9.1 | 61 to 63 | 890.7 | 4.2 to 4.7 from 249 (west southwest) | 9.4 to 9.7 | 37 to 45 | 1.101 |

Today's air was 4 % denser than September 19, which is worth about 2 % thrust in the aircraft's favour. The cold pack and the wind outweighed it. September 6 was also windy and that wingless flight has the worst yaw wander of the September set (±25°/s), without saturation because the pack was warm and the airframe light.

## 6. Battery evidence in detail

**Idle drop, like for like.** All four motors at 1184 µs with the same props:

| Log | Rest V | Idle V | Drop |
|---|---|---|---|
| 16 (Sep 19, 22 °C) | 48.03 | 47.77 | 0.26 |
| 17 (Sep 19, 22 °C) | 47.40 | 47.18 | 0.22 |
| 33 (Oct 18, 9 °C) | 48.82 | 47.96 | 0.86 |

**Hover sag.** September 19 flight 1 with wings: 47.4 V rest to 44.4 to 44.7 V at a 1550 to 1600 µs collective. October 18: 48.8 V rest to 40.3 to 40.6 V at a 1740 µs collective. Even allowing 40 % more current today, the per amp drop is more than double.

**Failsafe thresholds over time.**

| Parameter | Log 14 | Log 16 | Log 17 | Log 33 | Per cell at log 33 |
|---|---|---|---|---|---|
| BATT_ARM_VOLT | 46.8 | 46.8 | 44.0 (set at 1303 s) | 44.4 | 3.70 |
| BATT_LOW_VOLT | 44.4 | 44.4 | 44.4 | 42.0 | 3.50 |
| BATT_CRT_VOLT | 42.0 | 42.0 | 42.0 | 40.8 | 3.40 |
| BATT_FS_LOW_ACT | 0 | 0 | 0 | 0 | none |
| BATT_FS_CRT_ACT | 4 QLand | 4 | 0 (set at 1207 s) | 4 QLand | |
| BATT_FS_VOLTSRC | 0 raw | 0 | 0 | 0 | |
| BATT_LOW_TIMER | 10 s | 10 | 10 | 10 | |

**Current sensing.** BATT_MONITOR 4, BATT_VOLT_PIN 8, BATT_CURR_PIN 4, BATT_VOLT_MULT 18.18, BATT_AMP_PERVLT 36.36 (PM02 V3 defaults on POWER1). BATT2_MONITOR 0. Measured current 0.0 to 0.6 A in every log, so the shunt sees only the FC branch. This matches the SPH-E-003 "battery interface build" where the PM02 hangs off a 16 AWG pigtail and its output is unused. AB-08 (FCHUB V/I to POWER2) is the fix and has been open since October 16.

**KDE ESC voltage cutoff.** The UAS HVC manual lists the cutoff as OFF by default and recommends OFF for multirotors. If it had been ON at 3.2 V per cell it would have halved motor RPM at 38.4 V, which the log does not show. Low priority check in the KDE Device Manager.

## 7. Thrust margin arithmetic

KDE 5215XF-220 with the 21.5 x 7.3 DP prop, from the KDE performance sheet (measured at 1,067 m and 20 °C, density close to the field). ArduPilot ceiling 1898 µs (Q_M_SPIN_MAX 0.95 of 1100 to 1940). Interpolated between the 10S and 12S columns on a (duty × voltage)² shape, so treat as ±10 %.

| Bus V | Max thrust per motor at 1898 µs, kgf | Four motors, kgf | µs for 5.0 kgf per motor | µs for 6.0 kgf per motor |
|---|---|---|---|---|
| 48.0 | 9.0 | 36.1 | 1594 | 1670 |
| 44.5 (September hover) | 8.4 | 33.4 | 1632 | 1714 |
| 40.5 (October 18 hover) | 7.4 | 29.8 | 1685 | 1775 |

Reading September's 1546 to 1615 µs collectives against this table puts the September aircraft at roughly 4 to 4.5 kgf per motor. Today's 1738 µs reads as about 5.5 kgf per motor, inflated by the yaw fight and the climb and descent. The team's as flown masses would settle this. Whatever the exact mass, a CCW pair that must run 100 to 140 µs above the collective to hold heading reaches the 1898 µs ceiling at a collective of about 1760 µs, which is where today's hover sat.

## 8. What the logs give AB-03 (output map read back)

From the log 33 parameters, for the SPH-E-006 table:

| Output | Function | Notes |
|---|---|---|
| 3 | 70 throttle (pusher) | I/O port. Min 1100, max 1940. Never commanded in log 33 |
| 5 to 8 | 33 to 36, motors 1 to 4 | I/O port. SERVO5 to 8 max 1940, trim 1100 on 5, 1500 on 6 to 8 |
| 9, 10 | 4 aileron, 4 aileron | FMU port. SERVO9 1165 to 1800 trim 1565, SERVO10 1100 to 1950 trim 1500 |
| 13, 14 | 80 and 79, ruddervators | Over CAN, CAN_D1_UC_SRV_BM 12288 (bits 12 and 13). SERVO14 reversed, trim 1400, max 1850. SERVO13 trim 1550 |
| CAN_D1_UC_ESC_BM | 0 | Correct for the KDE PWM set |
| Q_M_PWM_MIN / MAX | 1100 / 1940 | Matches the KDE range |
| Q_RC_SPEED | 400 Hz | Within the KDE 50 to 600 Hz window |
| ARSPD_TYPE / BUS / USE / AUTOCAL | 1 / 1 / 0 / 0 | MS4525 detected, calibrated at boot, not yet used for control |
| GPS1_TYPE / GPS2_TYPE | 2 (u-blox serial, F9P) / 9 (DroneCAN, Here4) | September flew the Here4 alone as GPS1 |
| RNGFND1_TYPE | 27 (Benewake TF03 driver) | TF02-Pro reporting good status, 0.1 to 17.7 m seen. Driver choice to confirm with Zeynep |
| SERIAL2 | CRSF (RP4TD), link 250 Hz | |

## 9. Recommended order of work before the next hover

1. Airframe inspection after the hard touchdown: gear, boom roots, motor mounts, prop tips (finding 9).
2. Motor alignment survey with an inclinometer, static and strapped down at hover thrust, plus prop handedness and rotation per position. Shim. This is AB-12 and it is the item that decides whether the aircraft can hover in wind at all (finding 2).
3. Close AB-08 so the next log has current. Calibrate against a clamp meter (finding 4).
4. Pack: temperature before arming, cell checker before and after, bench capacity and resistance test, charger termination check (findings 3 and 6).
5. Re-enable the motor emergency stop on RC8 or document what replaces it (finding 10).
6. Parameters: Q_WVANE_GAIN back to 0 until the split is under 50 µs. Leave the battery thresholds alone until current works, then move to the sag compensated source. Export the parameter file (findings 5 and 7).
7. Next hover: warm pack, nose into wind, wind under 3 m/s, QLOITER, 2 to 3 m, 60 s steady so Q_M_THST_HOVER relearns, then read the CCW minus CW split and the pack sag from the log before anything else. Transition attempts wait until a calm hover shows yaw output under 0.1 and the collective under 1650 µs.

## 10. Questions for the team

- As flown mass on September 19 (with wings) and October 18 (with the pusher kit)?
- Pack temperature at arming today, and where the pack was kept before the flight?
- Charger model and termination voltage per cell for the ProFuse 22 Ah pack. Has the pack ever been balance charged and discharged on the bench to measure delivered Ah?
- Has anyone checked motor rotation direction and prop handedness per position since the September 21 build, and were the motor mounts touched when the pusher kit went on?
- Was the pack at 48.8 V by design (storage charge) or is that what the charger calls full?

## Method and sources

- pymavlink 2.4.49 DFReader over the five .BIN files. Per armed segment statistics on BAT, RCOU, QTUN, MOTB, RATE, ATT, VIBE, POWR, GPS, RCIN, RFND, PIQ*, MSG, ERR and PARM, plus 0.5 s time series around the October 18 hover. Scripts in the session scratchpad, not committed.
- Weather: Open-Meteo archive API (ERA5) for September 6 and 19, Open-Meteo forecast API with past days for October 18, hourly, Europe/Istanbul time zone.
- KDE5215XF-220 performance sheet (XLSX rev 2) and product page: https://www.kdedirect.com/products/kde5215xf-220
- KDE UAS55HVC product page and UAS ESC instruction manual rev PR13 (voltage cutoff default OFF): https://www.kdedirect.com/products/kde-uas55hvc
- ArduPilot motor mixer behaviour on saturation (roll and pitch, then yaw with Q_M_YAW_HEADROOM, then collective) per AP_MotorsMatrix.
- Quiver arm alignment precedent: Quiver Initial Configuration Guide §17.6 and §17.7, July 23 and 24, 2026.
