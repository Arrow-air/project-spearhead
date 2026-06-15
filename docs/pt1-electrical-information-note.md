# Spearhead PT1 Electrical Architecture Information Note

Internal doc ID: SPH-E-002. Companion to SPH-E-001 (Spearhead Electrical Master). AIP-005 submission ID to be assigned on GBC submission.

# Status

`Valid`

`Revision History: None`

`Replacement Log: None`

`Reference: SPH-E-001 (Spearhead Electrical Master), SPH-E-003 (PT1 parts list), AIP-006 (requirements)`

# Project Description

Project Spearhead is an Arrow Air fixed-wing QuadPlane UAV at ~25 kg MTOW, a Quad-X VTOL with an IC pusher for cruise. This information note records the Prototype 1 (PT1) electrical architecture: what was considered and what PT1 will use. It is a review summary. Full schematics, the bill of materials, and the work packages live in the electrical master document (SPH-E-001). The order-facing parts list with priorities and lead flags is SPH-E-003 (`docs/pt1-parts-list.md`).

PT1 covers the electric VTOL phase and its supporting avionics: flight controller, power distribution, battery and health monitoring, ESC and motor wiring, GPS and navigation, tail control-surface conversion, servos, the short-range RC and telemetry link, heading LEDs, pre-flight diagnostics, and the safety disconnects. The IC engine, generator, and long-range (>200 km) telemetry are Phase 2 and later and are out of scope here.

# Bounty or Grant Proposal Document

This information note documents an internal engineering work package (PT1 electrical planning and layout), not a discrete bounty or grant. Requirements derive from the Project Spearhead proposal AIP-006 §9. 

# Methodology

PT1 follows one principle: stay off the shelf and minimize custom PCBs. Each subsystem was chosen by weighing the final product requirement against what PT1 actually needs to achieve VTOL hover and transition. Where a final design choice adds mass or complexity that PT1 does not need, PT1 takes the simpler path, and this note records what stays reserved for the final design. Decisions stem from Spearhead engineering calls through June 11, 2026 and independent research, with key component specs verified against vendor pages in June 2026.

# Results and Deliverables

PT1 electrical decisions:

| Area | Considered | PT1 choice | Why |
|---|---|---|---|
| Flight controller | Pixhawk 6C, Pixhawk 6X | **Pixhawk 6X Standard Set v2A** | Selected June 11: the 6X IO PWM pins drive the wing flaperons directly, saving a second CAN-PWM node and shelf space. Set includes the PM02D HV and an M9N GPS (bench/backup). Order from Holybro EU, the TR listing is ~3× and out of stock |
| FC power | dedicated BEC, PDB BEC, Holybro power module | PM02D HV on an XT60 pigtail at the PDB battery pads | Single supply, POWER2 unused for PT1 |
| Power distribution | custom PDB, bus bars, FCHUB-12S | Matek FCHUB-12S V2 | Off-the-shelf, 4x70A per channel, on-board BECs, a 440A current sensor, a built-in VBat divider, and a PINIO-switchable 12V pad. No custom PDB |
| Battery | smart (CAN) vs dumb, single 12S vs 2x 6S | ProFuse Super Nano 12S 16Ah 60C, XT150 socket (TR stock) | Smart pack added ~1 kg, not worth it. 60C clears the 15C floor by 4x. Monitoring via the FCHUB 440A sensor plus the V2 VBat divider |
| Battery split | XT90 or AS150U Y-splitter vs pad tap | No splitter: the FC power module taps the FCHUB battery pads | Same electrical node as a Y with one fewer connector pair in the 320A path. The pack's XT150 is the ground safing disconnect |
| Regulators | per-rail BEC mix, BEC12S-PRO trio, PM12S-3 pair | 2× Matek PM12S-3 (per unit: Vx 5.25/6/8V at 15A + fixed 5V/4A + fixed 12V/4A, TR stock) | Tail unit: ruddervator Vx rail + 5V for adapter and Here4. Forward unit: flaperon Vx rail. Vx switches to 8V for the final servos with no hardware change. Peripheral 5V stays on the FCHUB BEC |
| Avionics location | nose bay vs shelf | Avionics shelf above the forward battery, top access (fuselage v5) | Easier PT1 wiring and access. Shelf reserves a rectangular footprint for a future integration PCB |
| ESC signal | daisy-chain vs hub | CAN1 to the four ESCs via a CAN splitter/hub | Off-the-shelf fan-out, no custom board |
| Tail control | Here4-as-converter vs dedicated adapter (CAN-L4-PWM or CAN-L431) | Matek CAN-L4-PWM adapter, Here4 GPS-only | Avoids the single point of failure of one Here4 doing GPS and tail PWM. Stocked in Turkey, CAN-L431 is the documented alternate |
| Flaperon control | local wing CAN-PWM node vs FC IO PWM | Pixhawk 6X IO PWM outputs 1–2 | The 6X selection made the wing CAN node unnecessary, saving a board and shelf space. The ~1–1.5 m wing PWM runs follow the EMI rules in SPH-E-001 §4.6 |
| Tail power | CAN 5V run vs local regulation | Local HV-to-5V BEC (servos + adapter + Here4) | Avoids 5V drop over the 2.5 m run under servo load. Tail stays on the HV kill only |
| Servos | 5V test vs 6-8.4V HV | 5V servos for PT1 | Hover and control-surface motion only. Simplifies wiring and drops both 7.8V BECs. Final design 6-8.4V |
| GPS | single vs dual, second Here4 vs F9P-class | Dual GPS confirmed. Primary unit TBD: second Here4 (matched pair, GPS yaw) or F9P-class/salvage (compass yaw) | Redundant GPS (AIP-006 §9.4). Decision flagged in the parts list (D1) before the order ships |
| Altimeter | Ainstein US-D1 (radar) vs Benewake TF03-180 (LiDAR) | Ainstein US-D1 | Radar tolerates prop wash dust during VTOL landing. Specs corrected June 2026 (110g, 100 Hz, UART or CAN). A unit is believed on hand in Turkey, confirm before ordering |
| RC + telemetry link | SIYI HM30 (shared with Quiver), ExpressLRS + SiK 915 or RFD900x | TBD, deferred to D3 / OQ-04 | PT1 is VLOS ≤5 km, so a short-range MAVLink + RC link covers E-REQ-08. ExpressLRS CRSF can passthrough MAVLink at short range and drop the separate telemetry radio. The receiver and radio lines in SPH-E-003 wait on this call. Upgrade to RFD900x before the Phase 3 extended-range runs |
| Heading LEDs | shared peripheral 5V vs dedicated feed | WS2812B strips on their own 5V feed | E-REQ-11. The high transient draw stays off the peripheral 5V rail so it cannot disturb GPS, RC, or telemetry. Driven from a FC GPIO via the ArduPilot NTF subsystem |
| Battery and health monitoring | smart-pack BMS/CAN telemetry vs sensor-based | FCHUB-12S 440A sensor + ArduPilot coulomb counting, PM02D HV as BATT2, plug-in cell checker on the ground | The dumb pack has no per-cell data. Vehicle current from the FCHUB sensor, pack voltage from the V2 VBat divider or the PM02D, SoC from coulomb counting. ESC voltage, current, and temperature stream over DroneCAN. Satisfies E-REQ-03 and E-REQ-10 |
| HV/LV kill | contactor, MOSFET anti-spark switch, manual disconnect, software E-stop | TBD, with a June 11 criterion: the HV kill must be switchable from the RC link | §9.3 requires HV and LV kills. Candidate analysis in SPH-E-001 §4.3. Any reduced PT1 kill is a recorded deviation. The pack's plain XT150 (no anti-spark) adds weight to the inline hardware options |

Deliverables:
- Electrical master document (SPH-E-001, Rev 0.16): full architecture, schematics, cost baseline BOM, and work packages.
- This information note (SPH-E-002).
- PT1 parts list (SPH-E-003, `docs/pt1-parts-list.md`): order priorities, candidate vendors, lead flags, and the open decision register D1–D5.
- Power and CAN topology diagrams (SPH-E-001 §2) and the PT1 harness connection table H1–H21 (SPH-E-001 §4.6), covering which cable goes where, wire types, and connectors.

# Remarks

- Open items gating final detail (decision register in SPH-E-003): primary GPS unit (D1), HV kill implementation against the RC-switchable criterion (D2), RC and telemetry link (D3). D4 closed June 11 (Pixhawk 6X, only the EU-vs-TR-restock sourcing call remains) and D5 closed with the ProFuse pack. The FCHUB voltage sense question closed with the V2 built-in divider. Still pending: final servo torque figures from the hinge-moment study.
- PT1 is a test platform. The 5V servos, the dumb battery, and any simplified kill switch are PT1-specific. The final design reinstates 6-8.4V HV servos on an 8V rail (same PM12S-3 units, Vx switched to 8V) and the full independent HV and LV kills per AIP-006 §9.3.
- The battery interface is the pack's plain XT150 (the ground safing disconnect, not anti-spark, see SPH-E-001 §4.3). The FC power module taps the FCHUB battery pads, so no Y-splitter exists in the HV path.
- Pre-flight diagnostics (E-REQ-10) ride on ArduPilot built-in pre-arm checks (GPS lock, compass, RC failsafe, ESC arming, battery, CAN node health, altimeter sanity). No added hardware for PT1.
- Safety gate before first hover: the AMPX ESCs cut all four motors if CAN commands are absent for >2 s, so CAN bus health is a pre-hover check, not a nicety. Detail in SPH-E-001 §4.13.
- Airspeed: PT1 reuses the pitot from storage, transducer interface still open (OQ-06). A backup I2C DLVR sensor is carried in SPH-E-003 in case the stored unit has no usable transducer.
- Per AIP-005, on GBC submission this note is placed in its project folder as `ID - Title` with a zero-padded ID from the Spearhead submission counter, and any pictures go in a `picture/` subfolder.
