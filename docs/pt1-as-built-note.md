---
title: "Spearhead PT1 and PT1.5 Electrical As-Built Note"
sidebar_label: "PT1 As-Built Note"
sidebar_position: 7
---

# Spearhead PT1 and PT1.5 Electrical As-Built Note

**Doc ID:** SPH-E-006
**Revision:** 1, October 17, 2026
**Author:** errrks
**Status:** As-built record. Items marked **Open** or **Check** still need a measurement or a parameter read-back.

Companion to SPH-E-001 (electrical master, Rev 0.23), SPH-E-002 (PT1 decisions, Rev 7), SPH-E-003 (PT1 parts list), SPH-E-005 (PT1 wiring diagram), SPH-E-007 (KDE ESC placement research note, July 27, 2026) and SPH-E-008 (PT2 scope and reuse map).

## Purpose

SPH-E-001 describes the PT1 design as it stood in June 2026: MAD motors on DroneCAN ESCs, Kingmax servos on all four surfaces at 8 V and 8.4 V, Radiolink RC, an Ainstein radar altimeter, and a single 12S pack. The airframe that hovered on September 21 and passed control surface tests on October 15 was built differently in several places, mostly because the MAD shipment never cleared Turkish customs. This note records the as-built state from Erick's October 15 and 16 report so the design docs can stay as history and the team has one place to read what is actually on the aircraft.

PT1 is the hover airframe. PT1.5 is the same airframe with the pusher fitted for transition attempts. Items that exist only for PT1.5 are marked.

## As-built configuration

| Item | SPH-E-001 design | As built | Status |
|---|---|---|---|
| VTOL motors and ESCs | MAD 150KV + AMPX 80A on DroneCAN, `CAN_D1_UC_ESC_BM = 0x000F` | **KDE 5215XF-220** motors + **KDE UAS55HVC** ESCs on PWM (confirmed Oct 16). ESCs mounted in the fuselage bay, phase leads extended with 14 AWG along the booms (Option B of the placement note). VTOL props **21.5 x 7.3 in** (Erick, Oct 17), the KDE catalogue match for the 5215XF-220 is the **KDE-CF215-DP** carbon two-blade set, 25.8 g per blade, M4 screws. Five ESC control leads land on a PWM breakout wired to the 6C **I/O PWM OUT** port (MAIN, ArduPilot outputs 1 to 8) | Flown Sept 21 |
| ESC control lead rail | None (DroneCAN) | 5 V rail on the ESC breakout, fed from the FCHUB-12S V2 5 V BEC pads through a 3-pin servo plug on a spare breakout header. KDE control inputs are opto isolated and accept 3 to 12 V | Flown |
| PWM breakout boards | None (DroneCAN) | Two **Holybro PWM EXT Boards** (SKU 18118), one on each 6C PWM OUT port. **Pin 1 (red, VDD_SERVO) removed from both 10-pin cables on October 16**, so each + row is fed only by its own regulator and the 6C servo rail is unpowered | Fix applied, readback pending (AB-01) |
| Pusher (PT1.5) | Electric pusher on the main pack, design intent | KDE 5215XF-220 (same motor as the VTOL set, per Erick "all" motors, confirm) + KDE ESC, **20 x 12 wooden propeller**, on its own **ProFuse 12S 15 Ah 15C** solid-state pack (brand confirmed Oct 17, 2,090 g, 170×60×96 mm per the motorobit listing) with an **XT90-S** connector, pack mounted in the aft fuselage near the pusher. Control lead on the same I/O breakout as the VTOL ESCs. Aft pack placement feeds the PT1.5 CG check | Fitting |
| Wing control surfaces | Flaperons, Kingmax CLS3015S at 8.4 V from the robocombo buck, 6C IO PWM 1 and 2 | **Ailerons**, no flap function. Decided Oct 16: plain ailerons through PT1.5, no flap mix. **Kingmax CLS3015S** on both (confirmed Oct 16) at **8.0 V**, robocombo buck pot set to 8.0 V. Servo leads on a second PWM breakout wired to the 6C **FMU PWM OUT** port (AUX, ArduPilot outputs 9 to 14) | Tested Oct 15 |
| Tail servos | Kingmax CLS3015S | **DSSERVO DS3115MG** on both ruddervators (corrected Oct 16, an earlier answer in this session said Kingmax). 4.8 to 7.2 V, 17.3 kg·cm at 6.5 V, 15.0 kg·cm at 7.2 V, 2.0 A stall, metal gear, 2BB, 25T, 50 to 330 Hz. This is the 6 V substitute foreseen in the June 29 build note. **8 V would exceed its rating**, which is why the tail Vx is at 6 V | Tested Oct 15 |
| Tail servo rail | PM12S-3 Vx at 8 V | PM12S-3 Vx bridged to **6 V**, matched to the DS3115MG. Re-bridging to 8 V is not an option while the DS3115MG is fitted | Tested Oct 15 |
| Tail control path | Matek CAN-L4-PWM node 12, ruddervators over DroneCAN | As designed. Node and both ruddervators validated with the FC in the loop | Tested Oct 15 |
| CAN bus layout | 6C CAN1, passive splitter at the boom roots, AMPX ESC nodes, tail | **Splitter removed.** Trunk runs direct from the 6C CAN1 to the CAN-L4-PWM solder pads. No ESC nodes on the bus. Two devices only: node 12 and the Here4 | As built |
| Tail CAN connection | Here4 then node, H13 daisy | The Here4 hangs off the node's JST-GH 4-pin CAN connector, so the Here4 is the physical end of the bus | As built |
| CAN termination | Not specified | FC end only. **Measured 120 Ω across CAN_H and CAN_L at the 6C CAN1 pins (Erick, Oct 16)**, which is one terminator on the bus, the 6C's own. CAN-L4-PWM 120R jumper open (factory). Here4 `CAN_TERMINATE` at 0 (factory) | **Fix**, see F2 |
| Primary GPS | u-blox F9P, serial, shelf | F9P on the nose shelf, serial GPS1 | As designed |
| Secondary GPS | Here4, tail, CAN node 11 | As designed | As designed |
| Airspeed | Pitot from storage, transducer open (OQ-06) | **Holybro PX4 Airspeed V1.3** carrying a TE **MS4525DO** (MS4525DO-DS5AI001DP, 1 psi differential, I2C address 0x28, 5 V supply, about 3 mA, 4.1 g) on the 6C I2C port, pitot mast on the nose (Erick, Oct 17). Closes OQ-06, the SPH-E-003 line 19 backup transducer is not needed | As built, `ARSPD_TYPE` 1 and `ARSPD_BUS` 1 to confirm in AB-03 |
| RC link | Radiolink AT9S Pro + R9DS, SBUS into RC-IN (D3) | **ExpressLRS 2.4 GHz, RadioMaster RP4TD receiver on TELEM2 over CRSF.** Radiolink not fitted. Gives RC and bidirectional link telemetry on one serial port | As built |
| Telemetry radio | Holybro SiK V3 on TELEM1 | As designed, SiK V3 on TELEM1. No band overlap with the 2.4 GHz ELRS link | As built |
| Rangefinder | Ainstein US-D1 radar (UART or CAN) | **Benewake TF02-Pro LiDAR on TELEM3**, serial. Ainstein not fitted | As built, `RNGFND1_TYPE` 19 (Benewake TF02) to confirm |
| Avionics bay temperature sensors (July 10 plan) | Two I2C sensors through a splitter | **Not fitted** | Dropped for PT1 |
| LV kill switch | Switch on the peripheral 5 V rail | **Not fitted** | Dropped for PT1 |
| FC supply and monitoring | PM02 V3 on POWER1, FCHUB V/I on POWER2 | As designed. POWER2 vehicle V/I test scheduled October 16 | **Check**, test pending |
| Battery | 12S 22 Ah semi-solid, QS8-S | As designed, the 12S 22 Ah semi-solid pack is the flown VTOL pack (confirmed Oct 16). Brand is **ProFuse**, motorobit is the vendor, same family as the pusher pack. Pusher pack separate, see above | Flown |
| HV kill | None, QS8-S unplug | None. With the pusher pack fitted, safing the aircraft takes **two** unplugs | Note |

## Validation log

| Date | Test | Result |
|---|---|---|
| Sept 21, 2026 | Four hover flights, KDE propulsion, Thomas piloting | Safe. Yaw and PWM imbalance from spar flex is the open flight issue |
| Oct 15, 2026 | Ailerons and ruddervators from the FC, ruddervators over CAN-L4-PWM | All four surfaces move correctly. CAN to PWM path validated with the DS3115MG tail servos, Kingmax ailerons validated on the FMU port |
| Oct 15, 2026 | Tail PM12S-3 at 6 V Vx plus 5 V logic rail | Working |
| Oct 16, 2026 | CAN bus resistance at the 6C CAN1 pins, unpowered | 120 Ω, one terminator on the bus (F2) |
| Oct 16, 2026 | FCHUB current and voltage on POWER2 | Not yet run, still open (AB-08) |

## Findings and required checks

### F1. Servo rail tie through the 6C (fix applied October 16, readback pending)

**Status:** Erick removed the red wire (pin 1, VDD_SERVO) from both 10-pin breakout cables on October 16, before the continuity check was run. The rails are now separate by construction. What remains is the powered readback at the end of this section. The analysis and the check steps stay here as the record of why.

Both 6C PWM OUT connectors carry **VDD_SERVO** on pin 1 (red). PX4's Pixhawk 6C hardware page lists them as one input: "Servo input: VDD_SERVO pin of FMU PWM OUT and I/O PWM OUT (0 V to 42 V undamaged)". Holybro's PWM extension boards and the generic 10-pin breakouts pass pin 1 straight to the servo header + row and carry no regulator of their own.

If both breakouts are plugged in with their 10-pin cables as supplied, the ESC breakout rail (FCHUB 5 V BEC) and the aileron breakout rail (robocombo at 8.0 V) are connected together through the 6C. The 8 V buck then back-drives the FCHUB 5 V output. The KDE control leads do not care (3 to 12 V input), which is why nothing has visibly failed, but the FCHUB BEC is being held above its own setpoint and the 6C servo rail is carrying wing servo current it was never meant to.

Evidence strength: strong but not a schematic. Holybro has not published the 6C schematic. The PX4 hardware page gives one absolute maximum rating for "VDD_SERVO pin of FMU PWM OUT and I/O PWM OUT", the Holybro ports page uses the same signal name on pin 1 of both connectors, the PX4 FMUv6C board config speaks of "the Servo rail" in the singular, and an ArduPilot forum user reports the 6C rail pins "all tied together". Every Pixhawk FMU design so far shares the rail. Treat it as probable until measured.

**Continuity check (settles it, two minutes)**

1. Aircraft fully unpowered. Both the VTOL pack and the pusher pack unplugged. USB out.
2. Unplug both 10-pin cables from the 6C.
3. Meter on continuity. Probe pin 1 (red wire position) of the FMU PWM OUT socket and pin 1 of the I/O PWM OUT socket on the 6C itself.
4. A beep, or under a few ohms, means one net inside the 6C and the rails are tied once both cables are plugged in.
5. Confirmation under power, cables plugged back in, aircraft powered: measure the + row on the ESC breakout. 5 V means the rails are separate. 8 V means the robocombo is reaching it through the 6C.

**What is happening if they are tied**

Each breakout is passive. Its + row is one copper strip joined to pin 1 of the 10-pin cable and to the + pin of every 3-pin header. The FCHUB 5 V feeds the ESC breakout + row and the robocombo 8 V feeds the aileron breakout + row. If the two pin 1 wires meet inside the 6C, the two + rows are one node and two regulators with different setpoints are wired output to output. The 8 V buck wins. The FCHUB 5 V BEC output is held about 3 V above its setpoint, and the KDE control leads, the 6C servo rail and both breakouts all sit at 8 V. The KDE inputs tolerate that, which is why nothing has announced itself.

**Fix: lift pin 1 on both cables**

Pin 1 is the only power conductor in the 10-pin cable. Pins 2 to 9 are PWM signals and pin 10 is ground. With pin 1 removed, the cable carries signal and ground only. Each breakout + row is then fed only by the regulator wired to it, and the 6C servo rail floats unpowered. That is fine because the 6C does not draw FC power from the servo rail and the PWM signal level comes from the FC's own 3.3 V or 5 V output driver, not from VDD_SERVO. If a 6C servo rail backup supply is wanted later, feed it from one source on purpose.

Methods, best first:

1. **De-pin the crimp.** At the 6C end of the cable, lift the small retention latch over position 1 of the JST-GH housing with a fine blade or pin and withdraw the red wire's crimp. Fold the wire back along the loom and heat shrink it. Reversible. Mark the housing.
2. **Cut the red wire** mid-cable and heat shrink both ends. Quick, not reversible.
3. **Cut the breakout trace** from connector pin 1 to the + row, if reachable. Hardest to inspect later, so only if the cable cannot be modified.

Do both cables. Lifting only the ESC side leaves the 6C servo rail at 8 V from the aileron breakout, within rating but pointless. Lifting only the aileron side leaves the ESC breakout backfeeding 5 V into the 6C servo rail. Neither is harmful, both are untidy. Two lifted pins means each regulator owns exactly one load.

Ground is unaffected. Pin 10 stays. The FCHUB, PM02 V3, robocombo and 6C already share ground through the power module and the PDB.

After the fix, repeat step 5 above. ESC breakout + row at 5 V, aileron breakout + row at 8.0 V, 6C servo rail pin 1 at 0 V.

### F2. CAN bus terminated at one end only

With KDE ESCs on PWM the CAN bus is short on node count but about 2.5 m long: 6C CAN1, trunk down the boom, CAN-L4-PWM pads, Here4 on the node's GH connector. The splitter is out. Erick measured 120 Ω across CAN_H and CAN_L at the 6C CAN1 pins on October 16, which confirms exactly one terminator on the bus, the 6C's own. The tail end is open. It works on the bench and it is marginal in the airframe with five ESCs switching beside the trunk.

Fix: set `CAN_TERMINATE = 1` on the Here4 through the Mission Planner DroneCAN parameter editor and reboot. The Here4 is the physical end of the bus, so this is the correct spot, and it needs no soldering. Alternative: bridge the 120R jumper on the CAN-L4-PWM if the Here4 stub stays under about 30 cm.

Verify: bus unpowered, meter across CAN_H and CAN_L at the 6C CAN1 pins, same spot as the October 16 reading. Expect about 60 Ω after the change. 120 Ω means the Here4 parameter did not take or the node needs a reboot.

### F3. Pusher pack grounding, safing and monitoring (PT1.5)

The KDE control input is opto isolated, so the pusher pack does not need a power ground bond for the PWM to work. The control lead's 5 V and GND are the FC side of the isolator.

Bond it anyway: one 16 to 18 AWG wire from the pusher pack negative to the FCHUB negative. It carries nothing in normal operation. It stops a 12S pack floating against a grounded airframe next to five other ESCs, which is a static and insulation stress problem waiting for humid weather, and it gives any fault current a path that is not a 26 AWG signal ground. Confirm the opto isolation on the actual ESC with a continuity check between the control lead GND pin and the pack negative lead before relying on it.

Safing: the QS8-S unplug only de-energizes the VTOL pack. The pusher pack needs its own connector pulled. Add it to the pre-flight and post-flight checklist in the ArduPilot reference guide.

Monitoring: the FC has no view of the pusher pack. **Decided Oct 16: no pusher pack monitoring on PT1.5.** Accepted risk. Pusher endurance is managed by timer and the pack is checked on the ground between flights. A BATT3 voltage divider into a spare 6C ADC pin remains the cheapest upgrade if PT1.5 flies longer transitions.

### F4. Aileron rail at 8.0 V, keep it

The Kingmax CLS3015S is rated to 8.4 V maximum. Running an adjustable pot at 8.4 V leaves no margin for pot drift, load regulation and light load rise. The torque difference is about 33 kg·cm at 8 V against 35 kg·cm at 8.4 V, around six percent. 8.0 V stays. Threadlock the pot and mark the setting.

### F5. Ruddervators on DS3115MG at 6 V

The tail flies the DSSERVO DS3115MG, not the Kingmax. Vendor figures: 14.5 kg·cm at 4.8 V, 17.3 kg·cm at 6.5 V, 15.0 kg·cm at 7.2 V, so about 16 kg·cm at the 6 V rail. The Kingmax at 8 V would have given about 33 kg·cm. Decided Oct 16: **PT1.5 flies with the DS3115MG at 6 V.** Re-bridging the PM12S-3 to 8 V is not available while this servo is fitted, since its rating stops at 7.2 V.

Before PT1.5 transition attempts, Zeynep's hinge moment check (WP-E06) must use about 16 kg·cm as the available torque at the horn, not 25.7 or 33. If the margin is thin, the fix is a servo swap to the two spare Kingmax and a Vx move to 8 V, with a re-test of the tail path. The 320 mm servo lead and the 50 to 330 Hz input range both suit the CAN-L4-PWM tail install.

### F6. ArduPilot output map must follow the hardware

SPH-E-001 assumes outputs 1 to 4 are DroneCAN motors, 1 and 2 are physical flaperon PWM, and `CAN_D1_UC_SRV_BM = 0x000C` sends outputs 3 and 4 to the tail node. The as-built wiring breaks that assumption: the I/O port now carries five physical ESC signals and the ailerons sit on the FMU port. Zeynep reconfigured the output map for the KDE PWM set before the September 21 hover (confirmed Oct 16), so motors 1 to 4 are already mapped on the I/O port. The channel numbers are not recorded anywhere yet, and whether the ruddervator CAN bitmask and the DroneCAN ESC bitmask were updated at the same time is unknown. The parameter file needs to be read back and the result written into this table. Expected shape:

| Output group | Physical port | Function |
|---|---|---|
| Outputs 1 to 4 (or the four used IO channels) | I/O PWM OUT | `SERVOn_FUNCTION` 33 to 36, VTOL motors 1 to 4 |
| One IO channel | I/O PWM OUT | `SERVOn_FUNCTION` 70, throttle (pusher, PT1.5) |
| Two FMU channels | FMU PWM OUT | `SERVOn_FUNCTION` 4 (aileron) each, or 77 and 78 if separate left and right trim is wanted. No flap function (decided Oct 16) |
| Two unused output numbers, not physically wired | Over CAN | `SERVOn_FUNCTION` 79 and 80 (ruddervators), with `CAN_D1_UC_SRV_BM` set to those two outputs. The bitmask must not overlap the physical ESC outputs |
| `CAN_D1_UC_ESC_BM` | | 0. No DroneCAN ESCs on PT1 |

KDE UAS ESCs use fixed PWM endpoints rather than a calibration routine. Set `Q_M_PWM_MIN` and `Q_M_PWM_MAX` to the KDE published range (1100 to 1940 µs per KDE documentation, verify on the UAS55HVC page) and confirm the ESC update rate matches `Q_RC_SPEED`. The ESC label reads 600 Hz OPTO.

## Open items

| ID | Item | Owner | Needed by |
|---|---|---|---|
| AB-01 | Pin 1 lifted on both cables Oct 16. Remaining: powered readback, ESC breakout + row 5 V, aileron breakout + row 8.0 V (F1). Optional: the unpowered continuity check between the two 6C pin 1 sockets, for the record | Erick | Next power up |
| AB-02 | Set Here4 `CAN_TERMINATE = 1`, reboot, verify about 60 Ω at the 6C CAN1 pins unpowered. Baseline 120 Ω measured Oct 16 (F2) | Erick, Zeynep | Before next flight |
| AB-03 | Zeynep to export the parameter file from the flown configuration. Read back `SERVOn_FUNCTION` 1 to 14, `CAN_D1_UC_SRV_BM`, `CAN_D1_UC_ESC_BM`, `Q_M_PWM_MIN/MAX`, `Q_RC_SPEED`, `ARSPD_TYPE` (expect 1, MS4525DO), `ARSPD_BUS` (expect 1, external I2C), `ARSPD_USE`, `ARSPD_AUTOCAL`. Record the channel assignment here (F6) | Zeynep | Before next flight |
| AB-04 | ~~Pusher pack connector~~ XT90, closed Oct 16. Brand ProFuse, 12S 15 Ah 15C, 2,090 g, 170×60×96 mm, recorded in SPH-E-003 Oct 17. Remaining: ground bond wire, two-unplug safing procedure in the checklist, aft pack position (station) to Zeynep for CG, mass is 2,090 g from the listing (F3) | Erick | PT1.5 fitting |
| AB-05 | ~~Pusher pack voltage monitoring method~~ Closed Oct 16: none on PT1.5, accepted risk, timer plus ground check. Add the pusher timer to the flight test plan | Zeynep | PT1.5 first transition |
| AB-06 | Ruddervator hinge moment at transition speed against about 16 kg·cm (DS3115MG at 6 V), WP-E06. Swap to the spare Kingmax at 8 V if thin (F5) | Zeynep | Before PT1.5 transition |
| AB-11 | Read back `SERIAL2_PROTOCOL` (CRSF for the RP4TD), `SERIAL4_PROTOCOL` and `RNGFND1_*` for the TF02-Pro on TELEM3, `RC_PROTOCOLS` | Zeynep | With AB-03 |
| AB-07 | ~~Airspeed sensor chip and I2C address~~ Closed Oct 17: Holybro PX4 Airspeed V1.3, MS4525DO at 0x28. Remaining: cover the pitot at every power up so the startup zero in `ARSPD_OFFSET` is clean, and plan the ratio calibration (`ARSPD_AUTOCAL = 1` through a few circuits, then back to 0) before the transition attempt | Erick, Zeynep | Before PT1.5 transition |
| AB-10 | ~~VTOL propeller size and make on the KDE 5215XF-220 motors~~ Closed Oct 17: **21.5 x 7.3 in**, two blade. KDE lists this size for the 5215XF-220 as the KDE-CF215-DP (M4 production version). Spares: KDE back-ordered to year end, the 3 blade option was dropped Oct 2 on KDE's 12S failure warning | Erick | — |
| AB-12 | Yaw and PWM imbalance from the Sept 21 hover. Still open Oct 16, no fix applied. Suspected spar flex. Root cause and mitigation before PT1.5 | Zeynep, Alperen | Before PT1.5 transition |
| AB-08 | POWER2 vehicle V/I test result, `BATT_AMP_PERVLT` and `BATT_VOLT_MULT` as calibrated. Not run as of Oct 16 | Erick | Before next flight |
| AB-09 | ~~Breakout board model for both PWM ports~~ Closed Oct 16: two Holybro PWM EXT Boards (SKU 18118). ~~Still to add to SPH-E-003~~ Added to the SPH-E-003 on-hand table Oct 17. No action left | Erick | — |

## Impact on the design documents

- **SPH-E-001:** §4.2 and §4.3 (DroneCAN propulsion), §4.4 (servo rails), §4.6 harness table H3, H8, H9, H11, H15, §4.10, §4.11 and the ArduPilot bitmasks all describe the design, not PT1 as built. Rev 0.22 adds a pointer to this note rather than rewriting them, because the DroneCAN design returns with the MAD set in PT2.
- **PT1 wiring diagram:** the propulsion half still draws DroneCAN ESC nodes and the H11 CAN runs. It is correct for the tail and servo rails after the October 14 update and stale for the ESCs and the pusher pack. Redraw after AB-03 fixes the channel map.
- **SPH-E-003:** ~~add the 12S 15 Ah pusher pack and its connector (AB-04)~~ done Oct 17, ProFuse 12S 15 Ah on XT90 in the on-hand table. The KDE set with the 21.5 x 7.3 in props and the two Holybro PWM EXT Boards went into the on-hand table October 17 (AB-09 and AB-10 closed), and the pitot row and the line 19 backup transducer were updated the same day for the MS4525DO.
- **SPH-E-001 OQ-06 and §4 airspeed table:** closed October 17 with the MS4525DO as built, pointer added in Rev 0.23.
- **Doc ID collision, resolved Oct 16:** SPH-E-005 was used three times. The wiring diagram keeps SPH-E-005. The ESC placement research note is now **SPH-E-007**. The PT2 scope and reuse map is now **SPH-E-008**, and the PCB project READMEs under `engineering/electronics/pcbs/` were updated to cite it. This note is SPH-E-006.
- **SPH-E-002 D3:** the RC side (Radiolink AT9S Pro + R9DS on SBUS) is superseded as built by ExpressLRS 2.4 GHz with a RadioMaster RP4TD on TELEM2. The telemetry side (SiK V3 on TELEM1) stands. The Ainstein US-D1 altimeter selection (June 11) is superseded as built by the Benewake TF02-Pro on TELEM3.
- **Naming:** the wing surfaces are ailerons on PT1. The design docs keep "flaperon" as history. New PT1 and PT1.5 text should say aileron.

## Sources

- PX4 Pixhawk 6C hardware page, servo rail absolute maximum rating: https://docs.px4.io/main/en/flight_controller/pixhawk6c.html
- Holybro Pixhawk 6C ports, FMU PWM OUT and I/O PWM OUT pinouts: https://docs.holybro.com/autopilot/pixhawk-6c/pixhawk-6c-ports
- PX4 FMUv6C board config, "the Servo rail" singular: https://raw.githubusercontent.com/PX4/PX4-Autopilot/main/boards/px4/fmu-v6c/src/board_config.h
- ArduPilot forum, user report that the 6C rail pins are tied together: https://discuss.ardupilot.org/t/pixhawk-6c-not-driving-motors-with-remote-control/99228
- PX4 CAN wiring, termination at both ends: https://docs.px4.io/main/en/can/
- Matek CAN-L4 family 120R jumper wording (CAN-L4-RC page): https://www.mateksys.com/?portfolio=can-l4-rc
- AP_Periph parameter list, `CAN_TERMINATE`: https://ardupilot.org/dev/docs/AP_Periph-Parameters.html
- Holybro CAN termination software toggle: https://docs.holybro.com/can-termination-resistance-software-toggle
- KDE Direct electronics FAQ, opto isolation and control lead voltage: https://www.kdedirect.com/pages/support-center-faq-electronics
- DSSERVO DS3115MG specifications (vendor listing): https://rees52.com/products/dsservo-ds3115mg-15kg-digital-servo-high-torque-metal-gear-digital-servo-motor-180-270-360-degree-rotation-digital-servo-25t-metal-gear-for-arduino-1-8-1-10-1-12-rc-cars-toy-robot
- Holybro MS4525DO digital airspeed sensor, specifications and I2C address: https://holybro.com/products/digital-air-speed-sensor-ms4525do
- Holybro MS4525DO and MS5525DSO setup, `ARSPD_TYPE = 1` for the MS4525: https://docs.holybro.com/peripherals/ms4525do-and-ms5525dso-i2c-airspeed-sensor
- ArduPilot airspeed sensor setup, `ARSPD_BUS`, startup zero, `ARSPD_AUTOCAL`: https://ardupilot.org/plane/docs/airspeed.html
- ProFuse 12S 15 Ah 15C solid-state pack, weight, dimensions, XT90-S (motorobit listing): https://www.motorobit.com/444v-12s-15000mah-15c-solid-state-lipo-battery
- KDE-CF215-DP 21.5 x 7.3 in propeller set, recommended for the KDE5215XF-220: https://www.kdedirect.com/products/kde-cf215-dp

## Revision history

| Rev | Date | Changes |
|---|---|---|
| 1 | October 17, 2026 | Airspeed sensor identified as built: Holybro PX4 Airspeed V1.3 with the MS4525DO at I2C 0x28. AB-07 closed, AB-09 closed (PWM EXT Boards added to SPH-E-003), AB-10 closed (VTOL props 21.5 x 7.3 in), pusher pack recorded as ProFuse 12S 15 Ah on XT90, AB-03 read-back list extended with the airspeed parameters, OQ-06 closed in SPH-E-001 and SPH-E-002, SPH-E-003 pitot rows updated |
| 0 | October 16, 2026 | Initial as-built record from Erick's October 15 and 16 report |
