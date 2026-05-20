# Open Questions — Project Spearhead

Questions requiring team input before electrical layout and procurement can be finalized. Organized by owner.

---

## Alperen — Structural / Layout

**Boom length**
The AMPX 80A ESCs ship with 12AWG/800mm power leads. If the boom is longer than ~700mm, the leads need extending with an XT90S splice at the boom-fuselage junction. What is the current target boom length for the forward and aft booms?

**Square CF tube dimensions**
ESCs will be externally mounted on the booms. The AMPX 80A body is 79.6×36×23.5mm with M2×4 mounting holes and a finned heatsink that needs airflow. What are the outer dimensions of the square CF tube, and is a bracket solution already in mind?

**Motor shaft diameter**
The FLUXER PRO prop hub requires a Ø10mm center bore and a Ø20mm M3×4 bolt circle. Can you confirm the V8013 PRO output shaft is 10mm and that a matching prop adapter plate is included, or will we need to source one separately?

**Nose bay dimensions**
Approximate L×W×H of the nose avionics bay — needed to finalize the physical layout (Pixhawk 6C, GPS, RC receiver, BECs, telemetry radio).

**Battery bay dimensions**
Approximate L×W×H available for the battery. A 22Ah 12S LiPo is roughly 180×65×75mm and ~5.5kg. Needs to be confirmed and CG position locked before the fuselage frame is finalized.

---

## Zeynep — Flight Mechanics / ArduPilot

**Control surface torque requirements**
Final minimum torque figures (kg·cm at 7.4V) for ruddervator servos and flaperon servos from the control surface sizing study. These are needed to finalize servo selection (WP-E06).

**Here4 CAN-to-PWM test result**
What were the results of the bench test during Thomas's Ankara visit? Pass/fail per the success criteria in `docs/bounty-here4-tail-servo-trial2.md`? This gates the tail PCB design (WP-E03).

**Here4 breakout PWM rail topology (continuity test)**
On the Here4 breakout board, are the V+ pins of all 8 PWM outputs electrically common (shared rail), and are the GND pins likewise common? This determines BEC wiring for the servo power supply.
- If shared: one V+/GND feed to any output channel powers all 8 channels (standard Pixhawk-style power rail). One BEC connection suffices.
- If isolated: each channel needs its own BEC feed or a jumper between rails.

Test procedure (board unpowered, multimeter in continuity mode):
1. Probe V+ of output 1 against V+ of outputs 2 through 8. Expect continuity on all.
2. Repeat for GND of output 1 against GND of outputs 2 through 8.
3. Probe V+ against GND of a single output. Expect open circuit (no continuity).
4. Document result and photograph the rail trace if visible.

This needs to be answered before the tail PCB power input footprint is finalized (WP-E03).

**ArduPilot servo function assignments**
Confirming V-tail ruddervator assignments: SERVO_FUNCTION 79 (RuddervatorLeft) and 80 (RuddervatorRight), broadcast over DroneCAN via `CAN_D1_UC_SRV_BM`? Any other servo channel assignments to note before the harness is wired?

**GCS preference**
Mission Planner or QGroundControl for ground station? Affects telemetry radio and RC link selection.

---

## Team / General

**Nose GPS**
The April 6 call discussed a second Here4 in the nose cone for dual GPS and moving-baseline yaw. Has this been confirmed, or is a different GPS module going in the nose? Options under consideration:
- Second Here4 (same as tail unit)
- Different DroneCAN GPS (Here3, Holybro, etc.)
- Single GPS only (tail Here4)

**Telemetry and RC link**
Are we continuing with the SIYI HM30 from Quiver? If yes:
- Does Spearhead have a dedicated unit, or would it share with Quiver?
- Is video downlink needed in Phase 1, or just MAVLink + RC?
If not HM30, what transmitter will pilots use? RC receiver selection depends entirely on this answer.

**Smart battery unit**
Smart battery confirmed. Which specific unit? Need to verify it broadcasts standard DroneCAN BatteryInfo messages and is available at 12S ≥ 22Ah capacity.

**HV kill switch type**
Two options under consideration:
1. High-current contactor (TE EV200 or Gigavac GX11) + key switch — positive lockout, ~300–500g, needs 12V coil supply
2. Anderson PP75 pull-disconnect with key-lock — lighter, simpler, no coil supply

**Altimeter selection**
Ainstein US-D1 (radar, 50m range, 12Hz, UART — robust to dust and debris) or Benewake TF03-180 (LiDAR, 180m range, 100Hz, UART — lighter, higher update rate)?

**Pitot tube**
The original Spearhead pitot tube was found in Alperen's drawer. What model is it and what is the output interface (analog, I2C, or UART)? Needed to confirm whether an additional transducer is required.

**Wing panel connector**
Molex Mini-Fit Jr (same family as Quiver) confirmed for the detachable wing panel interface, or a different standard?

**Payload bay connector**
Is there an Arrow platform standard forming for the payload interface (power + CAN/Ethernet), or will this be Spearhead-specific?

---

## Phase 2 (lower urgency — for awareness)

**Electric starter**
Any response from Pilot RC on the DLE35RA auto-starter? If unavailable, Phase 2 initial tests will use pull-start. Separate 2S or 3S LiPo for starter supply, or tap from main battery?

**IC engine ignition battery**
Shared with the main 12S VTOL battery, or a dedicated smaller pack? Affects HV bus and harness design for Phase 2.

---

*Last updated: May 2026 — Erick*
*Tracked in electrical master doc: `docs/electrical-master.md` (SPH-E-001)*
