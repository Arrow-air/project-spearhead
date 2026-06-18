# Spearhead PT1 Electrical Architecture and Decisions Information Note

Internal doc ID: SPH-E-002. Companion to SPH-E-001 (Spearhead Electrical Master) and SPH-E-004 (PT1 Electrical Planning and Design Foundation).

# Status

`Valid`

`Revision History: Rev 1 (June 17, 2026) — FC to full Pixhawk 6C + PM02 V3, battery to 12S 22Ah semi-solid (motorobit TR), primary GPS to salvaged serial F9P, servos to Kingmax CLS3015S HV on PT1, wing rail to a robocombo buck (tail PM12S-3 kept), FCHUB V/I to 6C POWER2 (Setup B). Rev 2 (June 17) — no hardware HV kill (software E-stop + QS8-S unplug, §9.3 deviation), final battery specs (QS8-S, 3,709g, 190×78×126mm), avionics compartment ~180×250×90mm, consistency pass. Rev 3 (June 18) — merged the standalone decisions-and-open-questions record into this note: adds the decision register (D1–D5), the deviation and accepted-risk register, the open-question list (formerly docs/open-questions.md), and pending confirmations. Tracks SPH-E-001 Rev 0.18.`

`Replacement Log: As of Rev 3 this note is the canonical open-question list, superseding docs/open-questions.md. The decision register in SPH-E-003 "Priority 3" remains the order-gating view; this note carries the engineering record.`

`Reference: SPH-E-001 (Spearhead Electrical Master), SPH-E-003 (PT1 parts list), SPH-E-004 (planning foundation)`

# Project Description

Project Spearhead is an Arrow Air fixed-wing QuadPlane UAV at ~25 kg MTOW, a Quad-X VTOL with an IC pusher for cruise. This information note records the Prototype 1 (PT1) electrical architecture and its decision state: what was considered, what PT1 will use, what is closed, what deviates from requirements, and what remains open. It is a review summary and a running status record. Full schematics, the bill of materials, and the work packages live in the electrical master document (SPH-E-001). The reasoning behind each choice lives in the planning foundation (SPH-E-004). The order-facing parts list with priorities and lead flags is SPH-E-003 (`docs/pt1-parts-list.md`).

PT1 covers the electric VTOL phase and its supporting avionics: flight controller, power distribution, battery and health monitoring, ESC and motor wiring, GPS and navigation, tail control-surface conversion, servos, the short-range RC and telemetry link, heading LEDs, pre-flight diagnostics, and the safety disconnects. The IC engine, generator, and long-range (>200 km) telemetry are Phase 2 and later and are out of scope here.

# Bounty or Grant Proposal Document

This information note documents an internal engineering work package (PT1 electrical planning and layout), not a discrete bounty or grant. Requirements derive from the Project Spearhead proposal AIP-006 §9. 

# Methodology

PT1 follows one principle: stay off the shelf and minimize custom PCBs. Each subsystem was chosen by weighing the final product requirement against what PT1 actually needs to achieve VTOL hover and transition. Where a final design choice adds mass or complexity that PT1 does not need, PT1 takes the simpler path, and this note records what stays reserved for the final design. Decisions stem from Spearhead engineering calls through June 17, 2026 and independent research, with key component specs verified against vendor pages in June 2026. The reasoning chain behind each trade is documented in SPH-E-004.

Decision state is tracked in four registers under Results: the considered-versus-chosen summary, the decision register (D1–D5, the procurement-gating calls), the deviation and accepted-risk register (departures from an AIP-006 SHALL and accepted single points of failure), and the open-question list (items needing team input, by owner). A decision is closed when the part and its rationale are fixed and only confirmation actions may remain.

# Results and Deliverables

## PT1 electrical decisions (considered versus chosen)

| Area | Considered | PT1 choice | Why |
|---|---|---|---|
| Flight controller | Pixhawk 6C Mini, full 6C, 6X | **Pixhawk 6C + PM02 V3 (full size)** | Reclosed June 17 (supersedes the June 11 6X). The 6C IO PWM still drives the flaperons (wing CAN node stays dropped), the full 6C gives 5 UARTs (vs the Mini's 4) for the serial F9P, and its two analog POWER ports carry the FCHUB vehicle V/I. No GPS in this bundle. Sourced from rx-dynamic TR |
| FC power | dedicated BEC, PDB BEC, Holybro power module | PM02 V3 (analog) on an XT60 pigtail at the PDB battery pads (POWER1) | Single supply on POWER1; POWER2 carries the FCHUB analog V/I sense, not a backup brick |
| Power distribution | custom PDB, bus bars, FCHUB-12S | Matek FCHUB-12S V2 | Off-the-shelf, 4x70A per channel, on-board BECs, a 440A current sensor, a built-in VBat divider, and a PINIO-switchable 12V pad. No custom PDB |
| Battery | smart (CAN) vs dumb, single 12S vs 2x 6S, ProFuse 16Ah vs solid-state 22Ah | motorobit 12S 22Ah 15C solid-state LiPo, 3,709g, 190×78×126mm, QS8-S (TR) | Closed June 17 (supersedes the ProFuse 16Ah). Dumb pack (smart added ~1 kg). 15C = 330A continuous, clears the floor and covers the 320A peak with more sag than 60C. 3,709g is ~113g lighter than the ProFuse, so no CG penalty for +6Ah. QS8-S is anti-spark. Monitoring via the FCHUB 440A sensor + V2 VBat divider |
| Battery split | XT90 or AS150U Y-splitter vs pad tap | No splitter: the FC power module taps the FCHUB battery pads | Same electrical node as a Y with one fewer connector pair in the 320A path. The pack's QS8-S is the ground safing disconnect |
| Regulators | per-rail BEC mix, BEC12S-PRO trio, PM12S-3 pair, PM12S-3 + buck | 1× Matek PM12S-3 (tail) + 1× robocombo DC-DC buck (wing) | Reclosed June 17. Tail PM12S-3: ruddervator Vx @ 8V + fixed 5V for adapter and Here4 (one TVS module for both rails). Forward PM12S-3 dropped: the wing flaperon rail moves to a robocombo buck @ 8.4V, since 8.4V gets full Kingmax torque and the PM12S-3 Vx caps at 8V. Peripheral 5V stays on the FCHUB BEC |
| Avionics location | nose bay vs shelf | Avionics shelf above the forward battery, top access (fuselage v5) | Easier PT1 wiring and access. Shelf reserves a rectangular footprint for a future integration PCB |
| ESC signal | daisy-chain vs hub | CAN1 to the four ESCs via a CAN splitter/hub | Off-the-shelf fan-out, no custom board |
| Tail control | Here4-as-converter vs dedicated adapter (CAN-L4-PWM or CAN-L431) | Matek CAN-L4-PWM adapter, Here4 GPS-only | Avoids the single point of failure of one Here4 doing GPS and tail PWM. Stocked in Turkey, CAN-L431 is the documented alternate |
| Flaperon control | local wing CAN-PWM node vs FC IO PWM | Pixhawk 6C IO PWM outputs 1–2 | The full 6C keeps the same IO PWM as the 6X, so the wing CAN node stays unnecessary. The ~1–1.5 m wing PWM runs follow the EMI rules in SPH-E-001 §4.6 |
| Tail power | CAN 5V run vs local regulation | Local HV-to-5V BEC (servos + adapter + Here4) | Avoids 5V drop over the 2.5 m run under servo load. Tail stays on the HV kill only |
| Servos | 5V test vs 6-8.4V HV | Kingmax CLS3015S HV (8.4V) on PT1 directly | Reclosed June 17: PT1 runs the final-class HV servos, dropping the 5V test stage. 35 kg·cm @ 8.4V, one SKU for flaperons + ruddervators. Note +100g aft tail mass (CG) |
| GPS | single vs dual, second Here4 vs F9P, CAN vs serial | Dual confirmed. Primary = salvaged u-blox F9P on a serial UART, secondary = Here4 on CAN | D1/OQ-12 closed June 17. Redundant GPS (AIP-006 §9.4). Mixed serial+CAN → compass yaw, no moving baseline. Salvaged F9P + on-hand antennas, no purchase |
| Altimeter | Ainstein US-D1 (radar) vs Benewake TF03-180 (LiDAR) | Ainstein US-D1 | Radar tolerates prop wash dust during VTOL landing. Specs corrected June 2026 (110g, 100 Hz, UART or CAN). A unit is believed on hand in Turkey, confirm before ordering |
| RC + telemetry link | SIYI HM30 (shared with Quiver), ExpressLRS + SiK 915 or RFD900x | TBD, deferred to D3 / OQ-04 | PT1 is VLOS ≤5 km, so a short-range MAVLink + RC link covers E-REQ-08. ExpressLRS CRSF can passthrough MAVLink at short range and drop the separate telemetry radio. The receiver and radio lines in SPH-E-003 wait on this call. Upgrade to RFD900x before the Phase 3 extended-range runs |
| Heading LEDs | shared peripheral 5V vs dedicated feed | WS2812B strips on their own 5V feed | E-REQ-11. The high transient draw stays off the peripheral 5V rail so it cannot disturb GPS, RC, or telemetry. Driven from a FC GPIO via the ArduPilot NTF subsystem |
| Battery and health monitoring | smart-pack BMS/CAN telemetry vs sensor-based; FCHUB sense on POWER1 (single-port) vs POWER2 (two-port) | **Setup B (two-port):** FCHUB-12S 440A sensor on 6C POWER2 as BATT1, PM02 V3 on POWER1 as BATT2, plug-in cell checker on the ground | The dumb pack has no per-cell data. Vehicle current + SoC from the FCHUB sensor (BATT1) on POWER2, redundant voltage from the PM02 V3 (BATT2) on POWER1, instances swapped so the real current sensor is primary. POWER2 is sense-only: FCHUB Cur + 1/21 VBat + a dedicated low-current G reference, no 5V. The 6C's two analog ports make this clean with no cutting the PM02 harness. Pins and scales in SPH-E-001 §4.13. ESC voltage, current, and temperature stream over DroneCAN. Satisfies E-REQ-03 and E-REQ-10 |
| HV/LV kill | contactor, MOSFET anti-spark switch, manual disconnect, software E-stop | **LV kill only; no hardware HV kill (closed June 17).** HV path: software E-stop + QS8-S anti-spark unplug for ground safing | §9.3 requires HV and LV kills, so PT1 carries a **recorded §9.3 deviation** (E-REQ-02): software E-stop replaces the HV hardware kill. Rationale: no clean contactor placement for PT1 and the QS8-S anti-spark handles connection inrush. Candidate analysis kept in SPH-E-001 §4.3; full HV kill returns for the final product |

## Decision register

### Closed

| ID | Decision | Closed | Outcome | Supersedes | Detail |
|---|---|---|---|---|---|
| D1 | Primary GPS (OQ-12) | June 17 | Salvaged u-blox F9P on a serial UART (GPS1), compass yaw, on-hand antennas. Secondary stays a single Here4 on CAN. Mixed serial-F9P + CAN-Here4 means no moving-baseline GPS yaw | Dual-Here4 / matched-pair options | SPH-E-001 §4.8, SPH-E-004 §4 |
| D2 | HV kill implementation (OQ-05) | June 17 | **No hardware HV kill for PT1.** Software E-stop (ArduPilot motor emergency stop on an RC switch, AMPX 2 s CAN watchdog backstop) + QS8-S anti-spark unplug for ground safing. LV kill retained. Recorded §9.3 deviation (see register below) | EV200 contactor, MOSFET anti-spark switch, RC-switchable criterion | SPH-E-001 §4.3, SPH-E-004 §3 |
| D4 | Flight controller | June 17 | **Full-size Pixhawk 6C + analog PM02 V3**, sourced rx-dynamic TR. IO PWM drives flaperons (wing CAN node stays dropped), 5 UARTs fit the serial F9P, two analog POWER ports carry the FCHUB vehicle V/I | June 11 Pixhawk 6X plan, 6C Mini, PM07/PM06 | SPH-E-001 §4.7, SPH-E-004 §4 |
| D5 | Battery format | June 17 | **Single motorobit 12S 22 Ah 15C semi-solid LiPo**, 3,709 g, 190×78×126 mm, QS8-S anti-spark socket. ~113 g lighter than the ProFuse despite +6 Ah. Only price left to confirm | June 11 ProFuse 16 Ah 60C, smart pack, 2×6S split | SPH-E-001 §4.1, SPH-E-004 §3 |

Sub-questions closed alongside the above:
- **FCHUB voltage sense:** closed by the FCHUB-12S V2 built-in 1K:20K divider (no external divider needed).
- **Battery-to-FC interface:** no Y-splitter; the FC power module taps the FCHUB battery pads (one fewer connector pair in the 320 A path). AS150U/XT150 Y kept as fallback only.
- **Servo class:** PT1 runs final-class Kingmax CLS3015S HV servos directly; the 5V test stage is dropped.
- **Bench FC = build FC:** the bench unit is the same full 6C as the flight controller, so Setup B two-port monitoring validates directly.

### Open

| ID | Decision | Options on the table | Owner | Needed by |
|---|---|---|---|---|
| D3 | RC link + telemetry radio (OQ-04) | SIYI HM30 (shared with Quiver?) vs ExpressLRS + SiK 915 / RFD900x. ExpressLRS CRSF can passthrough MAVLink at short range and drop the separate telemetry radio. Receiver and radio lines in SPH-E-003 wait on this | Team call | Before harness finalization. Upgrade to RFD900x-class before Phase 3 extended-range runs (E-REQ-07) |

## Deviation and accepted-risk register

| Item | Type | Requirement | What PT1 does | Rationale | Reversal path |
|---|---|---|---|---|---|
| No hardware HV kill | §9.3 deviation | E-REQ-02 (independent HV and LV hardware kills) | LV hardware kill retained; HV path uses software E-stop + QS8-S anti-spark unplug | No clean contactor placement for PT1 without added mass; QS8-S anti-spark handles connection inrush; software E-stop costs nothing | Full independent hardware HV kill for the final product; contactor/MOSFET analysis preserved in SPH-E-001 §4.3 |
| FC power on a single supply | Accepted SPOF | — (design margin, not a SHALL) | Pixhawk runs from one PM02 V3 on POWER1; POWER2 is used for the FCHUB V/I sense, not a backup brick | Keeps V1 power wiring simple; the two-port sense is worth more on PT1 than FC power redundancy | Move FCHUB sense to a spare ADC and put a backup BEC on POWER2 before extended flights |
| Dumb LiPo pack, no per-cell data | Accepted limitation | E-REQ-03 / E-REQ-10 (satisfied by sensor path) | No BMS/CAN; vehicle V/I from the FCHUB 440 A sensor + VBat divider; plug-in cell checker for ground checks | Smart pack added ~1 kg PT1 does not need | Smart/CAN pack or instrumented BMS if per-cell flight telemetry is required |
| Tail control path has no separate LV kill | Intentional | E-REQ-02 (LV kill scope) | Tail adapter, Here4, and ruddervators are HV-derived; de-energize only at QS8-S unplug | Preserves ruddervator authority whenever armed; resolves the adapter-loses-power-on-LV-kill concern | N/A for PT1; revisit with the full HV kill design |

## Open questions by owner

### Alperen — structural / layout

- **Boom length.** AMPX 80A ESCs ship 12 AWG / 800 mm power leads. Boom run exceeds 800 mm (confirmed June 11), so each lead is extended with an XT90S splice at the boom-fuselage junction. Need forward and aft boom target lengths to set extension length.
- **Square CF tube dimensions.** ESCs mount externally on the booms. AMPX body 79.6×36×23.5 mm, M2×4 holes, finned heatsink needs airflow. Need tube outer dimensions and the motor/ESC bracket plan.
- **Motor shaft diameter.** FLUXER PRO hub needs a Ø10 mm bore and Ø20 mm M3×4 bolt circle. Confirm the V8013 PRO output shaft is 10 mm and a matching adapter plate is included.
- **Avionics shelf dimensions (v5).** Current estimate ~180×250 mm footprint × ~90 mm tall between the longerons (June 12). Confirm in CAD and verify the set fits (Pixhawk 6C, PM02 V3, FCHUB-12S, CAN hub, primary F9P on the cover, RC receiver, telemetry radio) plus the reserved integration-PCB footprint.
- **Battery bay dimensions (forward bay, v5).** Pack is 190×78×126 mm, 3,709 g. Confirm the bay takes the 190 mm length plus strap, lead exit, and QS8-S clearance. Lock CG before the frame is finalized.
- **Tail harness routing (v5).** CAN trunk and the 18 AWG HV spur route along a tail boom. The tail PM12S-3 HV feed breaks out at that boom's ESC-lead extension joint, so the tail-carrying boom sets which ESC lead is tapped. Which boom, and external or internal? Internal routing waits on the June 5 structural check (holes in loaded CF members need analysis).

### Zeynep — flight mechanics / ArduPilot

- **Control-surface torque / hinge-moment.** Servos selected (Kingmax CLS3015S, 35 kg·cm @ 8.4V). Hinge-moment study still wanted to validate margin and set endpoint travel and linkage against the 25T horn (WP-E06). The +100 g aft from the two ruddervator servos feeds CG/stability.
- **CAN-to-PWM adapter validation (Matek CAN-L4-PWM).** One unit at the tail (node 12, ruddervators), `CAN_D1_UC_SRV_BM = 0x000C` broadcasting servo outputs 3–4. Confirm it enumerates as a DroneCAN servo node, accepts SERVO_FUNCTION 79/80, and outputs correct PWM under ArduPilot 4.x. Gates WP-E03/E04.
- **CAN-PWM adapter PWM rail topology (continuity test).** Are the V+ pins of the output channels common (shared rail), and likewise the GND pins? Determines whether one 8V Vx feed powers all channels. Procedure: board unpowered, probe V+ of output 1 against the others (expect continuity), repeat for GND, probe V+ against GND of one output (expect open), document and photograph the trace. Needed before tail BEC-to-adapter wiring (WP-E03).
- **ArduPilot servo function assignments.** Confirm SERVO_FUNCTION 79 (RuddervatorLeft) and 80 (RuddervatorRight) over DroneCAN via `CAN_D1_UC_SRV_BM`. Any other channel assignments to note before wiring?
- **GCS preference.** Mission Planner leads (Alperen and Thomas, including on Mac); QGroundControl is the alternative. Affects telemetry radio and RC link selection (feeds D3).
- **ESC sync freewheeling vs IPE prop auto-center (bench).** Does AMPX 80A sync freewheeling interfere with the V8013 PRO magnetic/mechanical auto-center indexing on spin-down? Confirm on the bench during first motor runs before locking ESC config.

### Team / general

- **Telemetry and RC link (feeds D3).** Continuing with the SIYI HM30 from Quiver? If yes: dedicated unit or shared with Quiver, and is video downlink needed in Phase 1 or just MAVLink + RC? If not HM30, which transmitter will pilots use? RC receiver selection depends entirely on this.
- **Pitot tube (OQ-06).** Original Spearhead pitot found in Alperen's drawer. Model and output interface (analog, I2C, or UART)? Determines whether the backup transducer is needed.
- **Wing panel connector (OQ-09).** Molex Mini-Fit Jr (same family as Quiver) confirmed for the detachable wing interface, or a different standard?
- **Payload bay connector.** Is an Arrow platform standard forming for the payload interface (power + CAN/Ethernet), or is this Spearhead-specific? (Phase 4 relevance.)
- **CAN splitter/hub model and termination (SPH-E-003 line 8).** Which passive JST-GH splitter for the CAN1 fan-out, and does it carry a built-in 120Ω terminator? Sets whether a separate terminator plug is needed and whether the boom stubs stay short enough at 1 Mbps. Fallback if stubs are marginal: daisy-chain the ESCs along the booms.
- **Per-branch ESC fusing (WP-E02).** Add the optional 100A fast-blow fuses per ESC branch on the FCHUB for PT1, or run unfused for the first prototype? Trade is fault isolation versus added resistance and connectors in the HV path; the 4×110A burst rating already covers the 320A peak.
- **Altimeter on-hand check.** Ainstein US-D1 selected (closed June 11). Confirm the believed-on-hand Turkey unit and its interface variant (UART vs CAN) before any order. Optional for first flight per June 12; an older Feather LiDAR may be a fallback (confirm CAN vs UART).

### Phase 2 (lower urgency, for awareness)

- **Electric starter.** Any response from Pilot RC on the DLE35RA auto-starter? If unavailable, Phase 2 initial tests use pull-start. Separate 2S/3S LiPo for the starter, or tap the main battery?
- **IC engine ignition battery.** Shared with the main 12S VTOL battery or a dedicated smaller pack? Affects the Phase 2 HV bus and harness design.

## Pending confirmations on closed decisions

These do not reopen a decision; they are the verification actions that remain.

| Item | Against | Status |
|---|---|---|
| Battery price | D5 (motorobit pack) | Confirm price; part and dimensions fixed |
| Forward-bay CG placement | D5 / layout | CAD; lock CG before the frame is built |
| Avionics-compartment fit check (~180×250×90 mm) | D4 / layout | CAD; verify the full set plus integration-PCB footprint |
| US-D1 on-hand + interface variant | Altimeter | On-hand check in Turkey before any order |
| Servo travel / linkage vs Kingmax 25T horn | Servos | WP-E06 hinge-moment + geometry |
| Bench validation of Setup B and the tail adapter | D4 / tail | WP-E04 on the full 6C bench unit |

## Deliverables

- This information note (SPH-E-002): the considered-versus-chosen summary, the decision register, the deviation and accepted-risk register, the open-question list, and pending confirmations. Canonical open-question list, superseding `docs/open-questions.md`.
- Electrical master document (SPH-E-001, Rev 0.18): full architecture, schematics, cost baseline BOM, and work packages.
- Planning and design foundation (SPH-E-004): the reasoning chain behind each PT1 choice and its final-product deferral.
- PT1 parts list (SPH-E-003, `docs/pt1-parts-list.md`): order priorities, candidate vendors, lead flags, and the order-gating decision register view (D2 and D3 open at issue, since closed except D3).
- Power and CAN topology diagrams (SPH-E-001 §2) and the PT1 harness connection table H1–H21 (SPH-E-001 §4.6), covering which cable goes where, wire types, and connectors.

# Remarks

- **Updating this note.** On a new call: record the outcome and date in the decision register, add any new deviation or accepted risk, and close or add the relevant open question. Then propagate to the affected SPH-E-001 section and, if procurement is touched, to SPH-E-003. Reasoning changes land in SPH-E-004.
- Only the RC and telemetry link (D3) remains open. D1, D2, D4, D5 closed June 17; the remaining work on closed items is the confirmations above, not re-decisions. The FCHUB voltage sense question closed with the V2 built-in divider.
- The §9.3 HV-kill deviation is the one departure from an AIP-006 SHALL on PT1. It is bounded (software E-stop + QS8-S unplug), recorded, and has a stated reversal path for the final product. Keep it visible in any PT1 review.
- PT1 is a test platform. The dumb battery and any simplified kill switch are PT1-specific. PT1 now runs the final-class Kingmax CLS3015S HV servos directly (flaperons @ 8.4V, ruddervators @ 8V), so the earlier 5V test-servo stage is gone. The final design still adds the full independent HV and LV kills per AIP-006 §9.3.
- The battery interface is the motorobit pack's QS8-S socket (the ground safing disconnect, anti-spark, see SPH-E-001 §4.3). The FC power module taps the FCHUB battery pads, so no Y-splitter exists in the HV path. Vehicle V/I rides the FCHUB sensor into the 6C's second analog POWER port.
- Pre-flight diagnostics (E-REQ-10) ride on ArduPilot built-in pre-arm checks (GPS lock, compass, RC failsafe, ESC arming, battery, CAN node health, altimeter sanity). No added hardware for PT1.
- Safety gate before first hover: the AMPX ESCs cut all four motors if CAN commands are absent for >2 s, so CAN bus health is a pre-hover check, not a nicety. Detail in SPH-E-001 §4.13.
- Airspeed: PT1 reuses the pitot from storage, transducer interface still open (OQ-06). A backup I2C DLVR sensor is carried in SPH-E-003 in case the stored unit has no usable transducer.
- Per AIP-005, on GBC submission this note is placed in its project folder as `ID - Title` with a zero-padded ID from the Spearhead submission counter, and any pictures go in a `picture/` subfolder.
