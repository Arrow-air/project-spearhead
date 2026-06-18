# Spearhead PT1 Electrical Planning and Design Foundation

Internal doc ID: SPH-E-004. Companion to SPH-E-001 (Spearhead Electrical Master) and SPH-E-002 (PT1 Electrical Architecture Information Note).

# Status

`Valid`

`Revision History: Rev 1 (June 18, 2026) — first issue. Captures the planning and engineering rationale behind the PT1 electrical architecture as it stood at SPH-E-001 Rev 0.18 / SPH-E-002 Rev 2 / SPH-E-003 (June 17). Where SPH-E-002 records what was considered and chosen, this note records why, and what each choice reserves for the final product.`

`Replacement Log: None`

`Reference: SPH-E-001 (master, Rev 0.18), SPH-E-002 (PT1 architecture note), SPH-E-003 (PT1 parts list)`

# Project Description

Project Spearhead is an Arrow Air fixed-wing QuadPlane UAV at ~25 kg MTOW, a Quad-X VTOL with an IC pusher for cruise. This is the foundation note for the electrical design: the planning, the framing, and the engineering reasoning that produced the Prototype 1 (PT1) architecture.

It is a thinking document, not a parts list and not a decision register. SPH-E-002 is the review summary of what was considered versus what PT1 will use, and it also carries the decision register, the deviation record, and the open-question list. SPH-E-003 is the order-facing parts list. This note sits underneath both and explains the design philosophy, the trade logic, and how each PT1 call ties back to a requirement and forward to the final product.

PT1 covers the electric VTOL phase and its supporting avionics. The IC engine, generator, BVLOS link, and payload bay are Phase 2 and later. They are out of scope for the PT1 architecture but in scope for the framing here, because most PT1 simplifications are only defensible once you state what the final product will need instead.

# Bounty or Grant Proposal Document

This note documents internal engineering planning for the PT1 electrical work package, not a discrete bounty or grant. Requirements derive from the Project Spearhead proposal AIP-006 §9.

# Methodology

## Design philosophy

PT1 follows one governing principle: **stay off the shelf and minimize custom PCBs.** PT1 is a test platform whose job is to validate VTOL hover, control-surface authority, structural behavior, and detachable-wing testing. It is not the product. Every subsystem was sized by weighing the final-product requirement against what PT1 actually needs to fly those tests. Where a final design choice adds mass, cost, or integration risk that PT1 does not need to retire, PT1 takes the simpler path and this note states what stays reserved for the final design.

Three consequences follow from that principle and recur through every decision below:

1. **A commercial board beats a custom board for V1.** A Matek FCHUB-12S replaces a custom PDB. A Matek CAN-to-PWM adapter replaces a bespoke tail node. A stock buck replaces a designed regulator stage. The reserved item in each case is the future integration PCB, for which the avionics shelf already holds a rectangular footprint.
2. **PT1 may carry deliberate single points of failure** that the final product will not, as long as each one is recorded and bounded. FC power on a single module and no hardware HV kill are the two live examples.
3. **PT1 runs final-class parts where running a throwaway test part would teach us nothing.** The clearest case is servos: PT1 flies the final HV servo class directly rather than a 5V test stage, because the test stage would not validate the torque, travel, or current the real surfaces see.

## Decision framework

Each subsystem trade was run the same way:

1. Pull the governing requirement from AIP-006 §9 (the traceability table is SPH-E-001 §3). Until the formal charter in proposal §5 is published, AIP-006 text is authoritative.
2. State what the final product needs to satisfy that requirement.
3. Ask what PT1 needs to satisfy it for the VTOL test campaign only.
4. If the two differ, take the simpler PT1 path and record the deferral. If they do not differ, build the final-class solution now.
5. Verify the chosen part against its vendor page (specs were checked in June 2026) and against the sourcing and customs rule before it enters SPH-E-003.

Decisions stem from Spearhead engineering calls through June 17, 2026 and independent research. The architecture moved materially between the June 11 and June 17 calls (FC, battery, regulators, GPS, kill switch), and this note reflects the June 17 state.

# Results and Deliverables

The reasoning is organized by the engineering question each cluster of decisions answers, not by part. For the resulting choices in table form see SPH-E-002. For the cable-level result see SPH-E-001 §4.6.

## 1. Requirements foundation and PT1 scope

AIP-006 §9 sets fifteen electrical requirements (E-REQ-01 through E-REQ-15). Seven are Phase 1 and gate PT1: VTOL thrust-to-weight (E-REQ-01), HV and LV kills (E-REQ-02), battery SoC telemetry (E-REQ-03), redundant GPS (E-REQ-04), radar or laser altimeter (E-REQ-05), short-range RF link (E-REQ-08), pre-flight diagnostics (E-REQ-10), heading LEDs (E-REQ-11), link-loss recovery (E-REQ-12), and a Pixhawk-standard FC with quadplane logic (E-REQ-13). The rest (obstacle avoidance, >200 km telemetry, payload bay, IC ignition, generator) belong to later phases and are deliberately excluded from PT1 so the test platform stays simple enough to build and debug.

The single largest scoping move is that **PT1 has no IC propulsion electronics at all.** No generator means no recharge-in-cruise path, no engine means no ignition or starter electrics. That removes the entire Phase 2 power-management problem from PT1 and lets the battery be sized purely for a hover-and-transition test budget rather than for endurance.

## 2. The PT1-versus-final-product framework

Every simplification below is paired with the final-product item it defers. This is the through-line of the whole architecture.

| PT1 simplification | Why acceptable for PT1 | Reserved for the final product |
|---|---|---|
| Dumb LiPo pack, no BMS or CAN | Smart pack added ~1 kg for per-cell data PT1 does not need; vehicle V/I comes from the FCHUB sensor | Smart/CAN pack or instrumented BMS if per-cell health is required |
| No hardware HV kill (software E-stop + QS8-S unplug) | No clean contactor placement for PT1, QS8-S anti-spark handles inrush; recorded §9.3 deviation | Full independent hardware HV kill alongside the LV kill |
| Commercial FCHUB-12S, no custom PDB | Off-the-shelf board carries the HV distribution, BECs, and current sensor | Custom integration PCB (footprint reserved on the avionics shelf) |
| FC on a single PM02 V3 supply | One supply keeps V1 wiring simple; recorded single point of failure | Backup BEC on a freed POWER2 (or FCHUB sense moved to a spare ADC) |
| Stock CAN-to-PWM tail adapter | Avoids a bespoke tail node and the Here4-as-converter single point of failure | Integrated tail node if the integration PCB absorbs it |
| Short-range RF link only | PT1 is VLOS ≤5 km, short-range MAVLink + RC covers E-REQ-08 | RFD900x-class link for the Phase 3 >200 km requirement (E-REQ-07) |

The discipline is that none of these are silent. Each one lives in the §9.3 deviation register or the open-question/SPOF list in SPH-E-002, so the gap between PT1 and the certifiable final product is always visible.

## 3. Power architecture reasoning

**Battery.** The energy driver for a QuadPlane is not sustained hover, it is the hover-and-transition budget, because the aircraft flies wing-borne for most of a mission. At ~85 A average hover draw a 22 Ah pack gives roughly 15 minutes of hover before reserve, far beyond any PT1 test flight. The first instinct in May was to drop to 16 Ah to save ~1.7 kg, and the June 11 call picked a 16 Ah ProFuse on that logic. June 17 reversed it: the motorobit 12S 22 Ah semi-solid pack is 3,709 g, which is ~113 g lighter than the 16 Ah ProFuse despite +6 Ah, because the cell density offsets the capacity. Once the weight penalty disappeared, the May energy trade no longer applied and capacity went back to 22 Ah. The C-rating reasoning is separate: 15C gives 330 A continuous, which clears the ≥15C floor and covers even the 320 A ESC-capped absolute peak, with the caveat that a 15C semi-solid pack will sag more under that peak than the 60C ProFuse would have, so the first high-throttle hover is a verification point.

**Distribution.** The FCHUB-12S decision is the design philosophy in one part: an off-the-shelf board with 4×70 A channels, on-board BECs, a 440 A current sensor, a built-in VBat divider, and a PINIO-switchable 12V pad, replacing a custom PDB and bus bars. The 70 A per-channel continuous rating comfortably clears the 80 A ESC ceiling and the 4×110 A burst covers the 320 A peak. The V2 built-in divider also closed the FCHUB voltage-sense question outright, removing a wiring unknown.

**The battery-to-FC tap.** The reasoning here is about connector count in the high-current path, not about voltage. A Y-splitter and a tap at the FCHUB battery pads create the same electrical node, but the pad tap puts one fewer connector pair in the 320 A path. So the AS150U/XT150 Y is dropped to a fallback and the PM02 V3 feeds off a 16 AWG pigtail at the pads. The QS8-S socket on the pack does the ground-safing disconnect job, and because it is anti-spark it manages mating inrush on its own, which is what later made the inline hardware kill unnecessary for PT1.

**Regulation.** The servo-rail plan moved from a symmetric "2× PM12S-3" idea to an asymmetric "PM12S-3 (tail) + buck (wing)" split for one concrete reason: the Kingmax CLS3015S reaches full 35 kg·cm torque at 8.4V, but the PM12S-3 Vx output caps at 8V. The tail keeps the PM12S-3 because it needs both an 8V servo rail and a 5V logic rail for the adapter and Here4 in one TVS-protected module, and 8V still gives ~33 kg·cm, which is ample for the ruddervators. The wing flaperons move to an adjustable buck set to 8.4V to recover the last bit of torque. The trade-offs of the cheap buck (lock the pot, confirm true synchronous 15A rating, scope for overshoot before connecting a servo) are accepted because the part is $15 and the risk is bench-checkable.

**Kill switches and the §9.3 deviation.** This is the most significant recorded compromise. AIP-006 §9.3 SHALL requires independent hardware HV and LV kills. PT1 keeps the LV kill but takes a recorded §9.3 deviation on the HV side: software E-stop (ArduPilot motor emergency stop on an RC switch, with the AMPX 2 s CAN watchdog as a backstop) plus unplugging the QS8-S anti-spark connector for ground safing. The reasoning chain: there is no clean contactor placement for PT1 without added mass and components, the EV200 contactor candidate is ~450 g and ~€149, the MOSFET alternative fails short and has no flight heritage, and the QS8-S anti-spark already manages the mating inrush that was the original argument for an inline kill. The software E-stop layer costs nothing and gets configured regardless. The full hardware HV kill returns for the final product, and the contactor/MOSFET analysis is preserved in SPH-E-001 §4.3 so the final-product work does not restart from zero.

## 4. Avionics and control reasoning

**Flight controller.** E-REQ-13 sets a Pixhawk-standard FC with quadplane logic as the floor. The interesting reasoning is why the full 6C beat both the Mini and the June 11 6X choice. Two drivers: the full 6C has 5 UARTs versus the Mini's 4, which the serial F9P primary GPS needs, and it has two analog POWER ports, which is what allows the FCHUB vehicle V/I to feed POWER2 directly. The analog PM02 V3 flipped from a constraint to an advantage in the process: the 6X's digital PM02D rode V/I over I2C and could not have fed an analog sense port, so the analog module is precisely what makes the two-port monitoring scheme work. The 6C IO PWM pins also drive the flaperons directly, which is what let the wing CAN-PWM node stay dropped.

**Flaperon control.** The wing CAN-PWM node considered on June 11 was dropped the moment the FC kept usable IO PWM pins. The cost is a ~1 to 1.5 m PWM run per side from the FC to each wing root, which is exactly the run the EMI rules in SPH-E-001 §4.6 exist to handle. The benefit is one fewer CAN node and one fewer thing to enumerate and debug. For a test platform, removing a node is worth managing a known wiring discipline.

**Tail control.** The choice of a dedicated CAN-L4-PWM adapter over using the Here4 as the PWM converter is a single-point-of-failure argument. One Here4 doing both GPS and tail PWM couples a navigation failure to a control-surface failure. Splitting them keeps the Here4 as GPS-only and puts ruddervator PWM on a $50 adapter that is in stock in Turkey, with the CAN-L431 as a documented alternate.

**Tail power topology.** The tail runs its own HV-to-5V regulation rather than a 5V run from the nose, because a 2.5 m 5V run under servo load would drop unacceptably. A consequence is that the tail control path is HV-derived and stays powered until the pack is unplugged. This is intentional: it preserves ruddervator authority whenever the aircraft is armed, and it resolves the earlier worry about the adapter losing power on an LV kill.

**Servos.** PT1 runs the final-class Kingmax CLS3015S HV servos directly and drops the earlier 5V test stage. The reasoning is that a 5V test servo would not exercise the torque, travel, or current the real HV surfaces see, so the test stage would teach nothing while still costing a build cycle. One SKU covers flaperons (8.4V) and ruddervators (8V). The cost is +100 g of aft tail mass, which feeds straight into the CG and stability work as a known input.

**GPS.** E-REQ-04 requires redundant GPS. The chosen pair is a salvaged u-blox F9P on a serial UART (primary) plus a Here4 on CAN (secondary). The serial-plus-CAN split is deliberate interface diversity, and the salvaged F9P with on-hand antennas costs nothing. The known limitation is that a mixed serial-plus-CAN pair cannot do moving-baseline GPS yaw, which needs two matched RTK units on one interface, so Phase 1 yaw comes from the compasses. The ~2 m nose-to-tail separation still gives spatial diversity for the redundant fix. If GPS yaw is wanted later, the path is two matched units on one interface.

## 5. Sensing and monitoring reasoning

**Health monitoring (Setup B).** The dumb pack has no per-cell data, so vehicle current and SoC have to come from a sensor rather than the battery. Setup B puts the FCHUB 440 A sensor on POWER2 as the primary monitor (BATT1) and the PM02 V3 on POWER1 as redundant voltage (BATT2), with the instances swapped so the real current sensor is primary. The reason this is clean rather than a hack is the 6C's second analog POWER port: it lets the FCHUB V/I feed a dedicated sense port with no cutting of the PM02 harness. ESC voltage, current, and temperature stream separately over DroneCAN. Together this satisfies E-REQ-03 and E-REQ-10. A plug-in cell checker covers per-cell ground checks that the dumb pack cannot report. The bench FC is the same full 6C as the build, so Setup B validates on the bench exactly as it runs.

**Altimeter.** E-REQ-05 allows radar or laser. The Ainstein US-D1 radar was chosen over the Benewake TF03 LiDAR because radar tolerates the prop-wash dust kicked up during VTOL landing, which is precisely the regime the altimeter matters most. A unit is believed on hand in Turkey, which would make it a zero-cost item, pending an on-hand and interface-variant check.

**Pre-flight diagnostics and the CAN watchdog.** E-REQ-10 is satisfied with ArduPilot built-in pre-arm checks (GPS lock, compass, RC failsafe, ESC arming, battery, CAN node health, altimeter sanity) and no added hardware. One item is not a nicety but a safety gate: the AMPX ESCs cut all four motors if CAN commands are absent for >2 s, so CAN bus health is a hard pre-hover check.

**Heading LEDs.** E-REQ-11. The only real design choice is to put the WS2812B strips on their own 5V feed rather than the peripheral rail, so their high transient draw cannot disturb GPS, RC, or telemetry.

**Airspeed.** PT1 reuses the pitot from storage, with the transducer interface still open. A backup I2C DLVR sensor is carried so an unusable stored transducer does not block the build.

## 6. Cross-cutting engineering principles

These shaped multiple decisions rather than one.

- **CG is the master constraint on placement.** The battery is the largest single CG item at 3,709 g and drives the forward bay (the 190 mm length is the binding dimension). The +100 g aft from the HV servos and the avionics-shelf-over-battery layout from fuselage v5 all feed the same CG accounting. Hardware mass is never chosen in isolation from where it sits.
- **Discrete wiring does not get Quiver's ground planes.** Quiver rejects this class of EMI mostly through multi-layer PCB routing. Spearhead's discrete harness does not have that, so separation, right-angle crossings, twisted pairs for PWM, and single-end shield grounding (at the FC) carry the load instead. This is why the flaperon PWM runs and the tail harness get explicit routing rules in SPH-E-001 §4.6.
- **Thermal margin is checked at the connector and the wire, not just the part.** The 8 AWG battery lead is short enough that the QS8 connector (~120 A continuous class), not the wire, is the thermal limit, and the worst-case 124 A hover draw sits at the top of that class, so connector temperature is a hover-test watch item. ESC power leads are matched 12-to-12 at the boom splice because the factory 12 AWG is rated to the full 80 A and Spearhead never sustains more than ~31 A/motor.
- **Accepted single points of failure are recorded, not hidden.** FC power on one module and the absence of a hardware HV kill are both deliberate and both written down with their reversal path. The standard is that a reviewer can always find the gap.
- **Sourcing follows the customs rule.** Buy in Turkey where stock exists, then EU stock, to keep customs simple, the same logic as the MAD order with the ATR certificate. This is why TR-side parts are sometimes accepted at a price premium (the PM12S-3 at ~4× Matek direct) when the customs simplicity is worth it.

## Deliverables

This note is the planning and rationale layer of the PT1 electrical doc set:

- SPH-E-004 (this note): design philosophy, decision framework, and the reasoning behind each PT1 choice with its final-product deferral.
- SPH-E-001 (master, Rev 0.18): full architecture, schematics, harness connection table, cost baseline BOM, work packages.
- SPH-E-002: review summary of considered-versus-chosen per subsystem, plus the decision register, the deviation register, and the open-question record.
- SPH-E-003: order-facing parts list with priorities, vendors, and lead flags.

# Remarks

- This note is intentionally stable. It records why the architecture is shaped the way it is. When a decision changes, the change and its new rationale land in SPH-E-002 and the affected SPH-E-001 section first, and this note is revised only if the underlying design philosophy or a major trade actually moves.
- The framing throughout is PT1-first with full-product context. Read a deferral here as a commitment recorded against the final product, not as a feature dropped. The reserved items (full HV kill, FC power redundancy, smart pack, integration PCB, long-range link) are tracked, not forgotten.
- The architecture moved substantially between June 11 and June 17 (FC 6X to 6C, battery ProFuse 16 Ah to motorobit 22 Ah, regulators to PM12S-3 + buck, GPS to serial F9P primary, HV kill to software-only). This note reflects the June 17 state at SPH-E-001 Rev 0.18.
- Per AIP-005, on GBC submission this note is placed in its project folder as `ID - Title` with a zero-padded ID from the Spearhead submission counter, and any pictures go in a `picture/` subfolder.
