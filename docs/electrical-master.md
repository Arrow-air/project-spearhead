# Spearhead Electrical Master Document

**Doc ID:** SPH-E-001
**Revision:** 0.1 — Draft
**Author:** errrks
**Date:** May 2026
**Status:** Phase 1 specified. Phases 2–4 are architecture stubs pending detailed design.

---

## Table of Contents

1. [Purpose and Scope](#1-purpose-and-scope)
2. [System Architecture](#2-system-architecture)
3. [Requirements Traceability](#3-requirements-traceability)
4. [Phase 1: Electric VTOL (Detailed)](#4-phase-1-electric-vtol-detailed)
   - 4.1 Battery System
   - 4.2 Power Distribution
   - 4.3 Kill Switches
   - 4.4 Voltage Regulation
   - 4.5 VTOL Motors and ESCs
   - 4.6 HV Wiring Harness
   - 4.7 Flight Controller
   - 4.8 GPS and Navigation Sensors
   - 4.9 CAN Bus Architecture
   - 4.10 Tail Electronics
   - 4.11 Servo System
   - 4.12 Telemetry and Communications
   - 4.13 Health Monitoring
   - 4.14 Heading LEDs
   - 4.15 Component Locations
   - 4.16 Phase 1 BOM and Cost Estimate
5. [Phase 2: Hybrid Integration (Stubs)](#5-phase-2-hybrid-integration-stubs)
6. [Phase 3: BVLOS (Stubs)](#6-phase-3-bvlos-stubs)
7. [Phase 4: Payload Bay (Stubs)](#7-phase-4-payload-bay-stubs)
8. [Work Packages](#8-work-packages)
9. [Open Questions](#9-open-questions)
10. [Revision History](#10-revision-history)

---

## 1. Purpose and Scope

This document defines the electrical architecture for Project Spearhead across all four development phases. Phase 1 (electric VTOL validation) is specified to component level. Phases 2–4 contain architecture stubs and open questions to be resolved as the project progresses.

Scope: power system, propulsion electronics, avionics, sensor integration, communications, harness design, component layout, and per-phase work packages.

Out of scope: structural design, ArduPilot parameter tuning (see `docs/reference-guide/ardupilot-quadplane.md`), flight mechanics.

---

## 2. System Architecture

### 2.1 Power Topology

```
                     ┌──────────────────┐
                     │  12S LiPo Pack   │
                     └────────┬─────────┘
                              │  AS150 / XT150
                     [HV Main Kill Switch]
                              │
                   ┌──────────┴───────────┐
                   │      HV Power Bus    │  44.4V nominal
          ┌────────┼───────┬──────────────┤
         ESC1    ESC2     ESC3           ESC4
      (DroneCAN)(DroneCAN)(DroneCAN) (DroneCAN)
          │        │       │              │
          M1       M2      M3             M4
       Front-R  Rear-L  Front-L        Rear-R

                   HV Power Bus
                        │
          ┌─────────────┼──────────┬──────────────────────────┐
          │             │          │                          │
     [5V BEC ×2]  [7.8V BEC]  [HV spur → Tail PCB]  [12V BEC, Phase 4]
          │        (CG bay,    (local 7.8V reg for          │
    Avionics rail  flaperons)   ruddervators + 5V        Payload rail
          │             │       for Here4 tail)
    [LV Kill Switch]    │
          │          Wing panel
     ┌────┴──┬────┐  connectors
  Pixhawk  Here4  Sensors/RC/
    6C    (nose)  Telemetry
```

### 2.2 Communications Topology

```
Pixhawk 6C
│
├── CAN1 ─── [120Ω term] ─── ESC1 ── ESC2 ── ESC3 ── ESC4
│                        ─── Here4 #1 (nose, node 10)
│                        ─── Here4 #2 (tail, node 11) ─── [120Ω term]
│
├── CAN2 ─── [reserved / Phase 3 expansion]
│
├── UART1 (TELEM1) ─── Ground telemetry radio (MAVLink 2.0)
├── UART2 (TELEM2) ─── [reserved / Phase 3 long-range link]
├── UART3 (GPS1)   ─── [unused, GPS on CAN]
├── UART4          ─── Radar/LiDAR altimeter
├── UART5          ─── RC receiver (CRSF or SBUS)
│
├── I2C   ─── External magnetometer (via Here4 internal, appears on CAN)
├── ADC   ─── Pitot static pressure transducer
│
├── PWM OUT (IO, outputs 1–8):
│     Output 1: Aileron (right flaperon)
│     Output 2: Aileron (left flaperon)
│     Output 3: RuddervatorLeft  (SERVO_FUNCTION 79, via DroneCAN to tail PCB)
│     Output 4: RuddervatorRight (SERVO_FUNCTION 80, via DroneCAN to tail PCB)
│     Outputs 5–8: VTOL motor commands via DroneCAN (not physical PWM)
│
└── USB-C ─── Ground configuration / parameter upload
```

VTOL motor commands go over DroneCAN via ESC_BM bitmask. Physical PWM outputs 5–8 are unused in Phase 1 unless DroneCAN falls back to PWM.

CAN bus length nose to tail: approximately 2.5 m. No signal integrity issue at 1 Mbps for this length.

### 2.3 Physical Location Overview

| Zone | Contents |
|---|---|
| Nose bay | Pixhawk 6C, Here4 #1 GPS, RC receiver, telemetry radio, 5V BECs (×2 redundant) |
| CG bay (under wing) | 12S smart battery (CAN telemetry), HV power bus / distribution, HV kill switch, 7.8V BEC for flaperon servos |
| Boom roots / tips (×4) | ESC externally mounted on boom near motor (within 150mm of motor centerline) |
| Motors (×4) | MAD V8013 PRO IPE 150KV, 3–5° forward tilt relative to aircraft body (motors vertical relative to gravity when aircraft is in cruise attitude) |
| Tail bay | Here4 #2 GPS, tail PCB (CAN-to-PWM + HV→7.8V regulation + 5V for Here4), ruddervator servo connectors |
| Wing panels (×2) | Flaperon servos, wing panel connector (signal + 7.8V + GND — 7.8V sourced from CG bay BEC) |
| Airframe exterior | WS2812 heading LED strips (wing tips, nose, tail) |

---

## 3. Requirements Traceability

Derived from the Project Spearhead DAO proposal (AIP-006, March 2026) and meeting decisions through May 2026. The formal project charter referenced in proposal §5 has not yet been published to the repo. Until then, the AIP-006 proposal text is the authoritative requirements source.

| ID | Requirement | Phase | Source |
|---|---|---|---|
| E-REQ-01 | VTOL thrust-to-weight ≥ 1.5:1 at MTOW 25 kg | 1 | AIP-006 §9.2 |
| E-REQ-02 | Independent hardware kill switches for HV and LV | 1 | AIP-006 §9.3 |
| E-REQ-03 | Real-time telemetry: battery SoC, generator output (Phase 2+), fuel level (Phase 2+) | 1/2 | AIP-006 §9.3, §9.5 |
| E-REQ-04 | Redundant GPS | 1 | AIP-006 §9.4 |
| E-REQ-05 | Radar or laser altimeter | 1 | AIP-006 §9.4 |
| E-REQ-06 | Front-facing obstacle avoidance sensor | 3 | AIP-006 §9.4 |
| E-REQ-07 | Long-range telemetry > 200 km | 3 | AIP-006 §9.5, §9.6 |
| E-REQ-08 | Secondary short-range RF link for local ops | 1 | AIP-006 §9.5 |
| E-REQ-09 | Payload bay: 12V power, CAN or Ethernet, quick-release | 4 | AIP-006 §9.7 |
| E-REQ-10 | Pre-flight diagnostics: battery, GPS, altimeter, sensor health | 1 | AIP-006 §9.8 |
| E-REQ-11 | Heading indicator LEDs | 1 | AIP-006 §9.8 |
| E-REQ-12 | Link-loss and engine-out autonomous recovery | 1/2 | AIP-006 §9.8, §2.5 |
| E-REQ-13 | Pixhawk-standard FC with quadplane transition logic | 1 | AIP-006 §9.4 |
| E-REQ-14 | IC engine ignition and optional starter | 2 | AIP-006 §9.2 |
| E-REQ-15 | Generator: avionics power and battery recharge during cruise | 2 | AIP-006 §9.3 |

---

## 4. Phase 1: Electric VTOL (Detailed)

Phase 1 scope: all systems needed for VTOL hover, control surface authority, structural validation, and detachable wing testing. No IC engine. No generator.

---

### 4.1 Battery System

#### Sizing

Based on the MAD V8013 PRO 150KV datasheet (measured at 48V, 25°C, sea level with FLUXER PRO 26×7.8 MATT prop):

**Datasheet voltage vs. 12S operating voltage:**

A 12S LiPo charges to 50.4V (4.2V/cell) and the datasheet was measured at 48V — placing the test conditions squarely within the normal mid-flight range of a 12S pack. The datasheet numbers are therefore directly representative rather than requiring large corrections. As the battery discharges, bus voltage drops and the motor must spin faster (more throttle) to maintain the same thrust, increasing current draw. The motor's operating zone shifts progressively across a flight.

Required thrust per motor: 25,000 g / 4 = **6,250 gf**

**Hover operating point across battery discharge (sea level):**

| Battery state | Bus voltage | Throttle (approx.) | Current/motor | Total draw | Motor zone |
|---|---|---|---|---|---|
| Full charge | 50.4V | ~57% | ~17A | ~68A | Continuous |
| ~75% SoC (matches datasheet) | 48V | ~60% | ~19.8A | ~79A | Continuous |
| ~50% SoC | 45V | ~63% | ~22A | ~88A | Continuous |
| Nominal | 44.4V | ~65% | ~24.4A | ~98A | Continuous / edge |
| Low battery warning | 42V | ~68% | ~28A | ~112A | Short-term zone |

Motor operating zones from datasheet footnote (motor thermal limits, not ESC limits — see §4.13):
- **Below 24A/motor:** Sustainable continuous
- **24–78A/motor:** Short-term only, ~10–30s
- **Above 78A/motor:** Non-working

The motor operates in the continuous zone for most of the usable flight window. It only enters the short-term zone near the low-battery warning threshold (42V / 3.5V/cell) — which is also the point at which landing should already be initiated. At altitude (1,200m, ~12% lower air density), add approximately 2–3% throttle to each row above, nudging the nominal and low-battery points further into short-term territory.

**Battery capacity sizing** (design case: 10-minute hover at average ~85A across full→low charge range):

- 10-minute hover: 85A × 10/60 = 14.2 Ah
- 20% reserve: 14.2 / 0.8 = **17.7 Ah minimum → 22,000 mAh recommended** (provides margin for altitude, repeated hover cycles, and transition phases)

**Peak current** (safety-margin thrust at 75% throttle, 36.6A/motor × 4):
- 4 × 36.6A = **146A** (in short-term zone, < 30s)
- At 22,000 mAh: 146 / 22 = **6.6C** (well within standard LiPo ratings)

**Thrust-to-weight at max thrust** (100% throttle, 15,032 gf/motor × 4 = 60,128 gf = 61.3 kg):
- T/W = 61.3 / 25 = **2.45:1** (E-REQ-01 of 1.5:1 exceeded)

#### Specification

| Parameter | Target |
|---|---|
| Chemistry | LiPo, high-discharge |
| Cell count | 12S |
| Nominal voltage | 44.4V |
| Full charge | 50.4V |
| Low-battery cutoff | 3.5V/cell → 42V pack (ArduPilot warning), 3.3V/cell → 39.6V (critical) |
| Capacity | **22,000 mAh** (minimum; 22,000 recommended) |
| C-rating | ≥ 15C continuous (provides 330A from 22Ah pack — ample margin for 146A safety-thrust draw) |
| Connector | AS150 (preferred) or XT150 |
| Form factor | Single pack, placed at CG — lock position with Alperen before fuselage is built |
| Weight (est.) | 5–6.5 kg depending on capacity and brand — largest single CG impact item |

**Smart battery confirmed.** A DroneCAN smart battery will be used. Specific unit TBD (WP-E01). Implications:
- External current sensor and power module (PM07) are **not required** — the battery broadcasts voltage, current, SoC, and temperature over CAN
- No balance tap wiring — cell balancing is handled internally by the BMS
- ArduPilot parameter: `BATT_MONITOR = 8` (UAVCAN/DroneCAN BMS)
- Battery connects to the HV bus via main power leads only (no balance harness)
- Verify the selected unit supports DroneCAN BatteryInfo messages and is available in a 12S form factor at ≥ 22Ah capacity

---

### 4.2 Power Distribution

| Parameter | Value |
|---|---|
| Architecture | Star topology: battery to central bus, individual runs to each ESC |
| Bus | Copper bus bars or high-current PCB (PCB preferred for Prototype 1 per Erick) |
| Bus current rating | 300A continuous (safety margin above 320A peak fuse protection) |
| Per-branch protection | 100A fast-blow fuse per ESC branch (optional for Prototype 1; assess during WP-E02) |
| Voltage monitoring | Via smart battery DroneCAN telemetry — no external power module required |

Candidate power distribution boards: Matek FCHUB-12S (12S, screw terminals, 200A rated) or custom PCB. Holybro PM07 is not needed since voltage and current monitoring come from the smart battery over CAN.

---

### 4.3 Kill Switches

| Switch | Covers | Does NOT cover | Placement | Current rating |
|---|---|---|---|---|
| HV Kill | 12S battery → ESC bus, tail PCB HV spur | — | Exterior, reachable within 3 s | ≥ 300A peak capability |
| LV Kill | 5V avionics rail (nose bay), 7.8V flaperon BEC (CG bay) | Tail PCB + ruddervator servos (on HV spur) | Exterior, separate from HV | ≥ 10A |

The tail PCB and ruddervator servos are powered directly from the HV bus spur and have no separate LV kill switch. They are de-energized only when the HV kill switch is thrown. This is intentional — ruddervator authority is required whenever the aircraft is powered.

HV kill switch options (decide in WP-E02):

| Option | Description | Pros | Cons |
|---|---|---|---|
| High-current contactor + key switch | TE EV200 100A contactor or Gigavac GX11 with external key | Positive lock, rated for repetitive switching | Heavy (300–500g), needs 12V coil supply |
| Anderson PP75 pull-disconnect with key-lock | Fused disconnect with mechanical latch | Light, simple, no coil supply | Requires extraction force; less positive lockout than contactor |
| E-Stop pushbutton + contactor | Mushroom-head button latches contactor open | Intuitive safety action | Contactor weight same as Option 1 |

LV kill switch: illuminated rocker or key switch, 10A @ 12V. Interrupts the nose 5V BEC output and the CG bay flaperon BEC only.

---

### 4.4 Voltage Regulation

| Rail | Voltage | Continuous Current | Loads | BEC type |
|---|---|---|---|---|
| Avionics | 5.0–5.2V | 8–10A | Pixhawk 6C, 2× Here4, RC receiver, telemetry radio, sensors, LEDs | Switching (efficiency critical at 12S input) |
| Avionics (backup) | 5.0–5.2V | 3–5A | Pixhawk 6C POWER2 input only | Second switching BEC or PM07 backup rail |
| Flaperon servo | 7.8V | 2–3A | 2× flaperon servos on wing panels | Switching BEC, located in CG bay. Short runs to wing panel connectors only — avoids long LV servo runs through fuselage. |
| Ruddervator servo | 7.8V | 2–3A | 2× ruddervator servos in tail bay | Regulated locally on tail PCB from HV bus spur. No 7.8V cable run through fuselage. |
| Payload | 12V | TBD | Phase 4 only — not needed for Phase 1 | — |

Redundancy: Pixhawk 6C accepts two independent power inputs (POWER1 and POWER2). Run primary BEC to POWER1 and a backup BEC to POWER2 for automatic hot-standby. Both source from the same HV bus but through separate regulators.

BEC candidates (decide in WP-E01): Matek MBEC6S (5–9V adjustable, 10A, 12S capable), Castle BEC Pro (20A, 6–12S), Holybro PM07 integral BEC (5.2V/6A, but rated only 60A on the HV side — need a separate high-current distribution board if used this way).

---

### 4.5 VTOL Motors and ESCs

**Status: ordered** (6 of each, MAD Motor Poland, expected delivery 2–4 weeks from May 1, 2026; ATR certificate requested for Turkish customs — zero import tax if provided).

#### Motors

| Parameter | Value |
|---|---|
| Model | MAD V8013 PRO IPE 150KV |
| KV | 150 |
| Weight | **591 g** |
| Battery | 12S (44.4V nominal; datasheet tested at 48V) |
| Max temperature | 87°C |
| Recommended hover thrust | 4.5–6.5 kgf (Spearhead at 6.25 kgf/motor — top of range) |
| Maximum thrust | **15 kgf** (with FLUXER PRO 26×7.8 MATT at 48V) |
| Electrical efficiency | > 80% (peak ~82% at 65–75% throttle) |
| Max propeller | 26" (boom spacing constraint) |
| Auto-center feature | Magnetic/mechanical indexing: props stop parallel to boom when power removed |
| Mount pattern | 40mm bolt circle, M4×4 |
| Bearings | Dual 6901 deep-groove |
| Quantity ordered | 6 (4 installed, 2 spares) |
| Motor tilt | 3–5° forward relative to aircraft body axis. The aircraft's natural cruise attitude (nose slightly elevated from wing incidence) places the motor thrust axes approximately vertical relative to gravity during hover — the structural tilt and aircraft attitude cancel. ArduPilot level calibration must be performed with the aircraft in cruise attitude (not flat on the ground) per `docs/reference-guide/ardupilot-quadplane.md §12`. |

**Performance table** (MAD datasheet, 48V, 25°C, sea level, FLUXER PRO 26×7.8 MATT):

| Throttle | Current (A) | Input Power (W) | Thrust (gf) | Efficiency (gf/W) | Elec. eff. (%) | Zone |
|---|---|---|---|---|---|---|
| 30% | 3.5 | 169 | 1,694 | 10.0 | 68.6 | Continuous |
| 40% | 6.8 | 326 | 2,894 | 8.9 | 75.1 | Continuous |
| 50% | 11.9 | 571 | 4,355 | 7.6 | 79.5 | Continuous |
| 55% | 15.6 | 745 | 5,313 | 7.1 | 81.0 | Continuous |
| 60% | **19.8** | **945** | **6,440** | **6.8** | 82.0 | **Continuous limit ≈ here** |
| 65% | 24.4 | 1,163 | 7,318 | 6.3 | 82.2 | Short-term (10–30s) |
| 70% | 30.0 | 1,420 | 8,476 | 6.0 | 82.3 | Short-term |
| 75% | 36.6 | 1,726 | 9,759 | 5.7 | 82.1 | Short-term |
| 80% | 41.7 | 1,961 | 10,389 | 5.3 | 82.0 | Short-term |
| 90% | 58.5 | 2,718 | 12,985 | 4.8 | 80.4 | Short-term |
| 100% | 77.8 | 3,569 | 15,032 | 4.2 | 78.0 | Short-term / non-working |

**Datasheet operating zones:**
- **Below 24A (< ~62% throttle at 48V):** Sustainable continuous
- **24–78A (~62–99% throttle):** Short-term only, ~10–30 s
- **Above 78A:** Non-working zone

At Spearhead's MTOW and nominal battery voltage (44.4V), hover falls near 65% throttle at sea level and ~69% at 1,200m altitude — placing each motor at or just above the continuous zone boundary. Monitor ESC temperatures during all hover tests. Set `Q_M_THST_HOVER` accurately to minimize steady-state throttle.

**Motor wire enamel risk:** If motor phase leads are shortened, inspect strand insulation for enamel or clear varnish coating. Remove at the termination point before soldering. A coated termination creates a high-resistance joint that will overheat under load. (Ref: Feather PC3 incident — see Erick's failure notes.)

#### ESCs

| Parameter | Value |
|---|---|
| Model | MAD AMPX 80A V2 (5–14S) DroneCAN |
| Protocol | DroneCAN V2.0 (also supports CyphalCAN; also accepts PWM 1050–1940µs) |
| Continuous current | 80A (requires good heat dissipation) |
| Current limit (burst) | 84A |
| Input voltage | 5–14S, 16–64V |
| BEC output | **5V / 200mA** — insufficient for avionics; use only as CAN signal reference or auxiliary. Dedicated BECs required for avionics and servo rails. |
| CAN telemetry fields | Voltage, current, temperature, operation status (broadcast in real time over DroneCAN) |
| PWM input voltage | 3.3V or 5V compatible |
| Protection level | IP67 (fully sealed, waterproof) |
| Operating temperature | -20°C to +65°C ambient |
| Thermal protection | >125°C: output reduces to 50% max. >140°C: full shutdown; will not restore until throttle zeroed AND ESC cools to <80°C |
| Short circuit protection | Shuts down output permanently — **requires power cycle to recover**. One motor loss per flight is unrecoverable in flight. |
| Throttle loss protection | Cuts power if CAN command absent for **>2 seconds**. ArduPilot CAN bus must maintain continuous ESC commands or motors cut. |
| Startup protection | Motor must start within 10s of throttle advance or ESC shuts down |
| Voltage protection | Below 16V or above 64V: alarm and no start (disabled during flight). 12S operational range 36–50.4V is within safe limits. |
| Sync freewheeling | Supported — active energy recovery on throttle reduction. Test interaction with IPE auto-center feature. |
| Weight | ~90g (excluding wires) |
| Dimensions | 79.6 × 36.0 × 23.5mm |
| Mounting | M2×4 screw holes (4 points) |
| Power wire | 12AWG / 800mm |
| Phase wire | 14AWG / 150mm |
| Signal wire | 1000mm (bare leads: Black=GND, White=PWM, Red=CANH, Green=CANL + 5V BEC wire) |
| Quantity ordered | 6 (4 installed, 2 spares) |
| Unit price | $60 (DroneCAN V2) + $30 shipping |

DroneCAN node ID assignment:

| Node ID | Motor position | Rotation | ArduPilot servo output |
|---|---|---|---|
| 1 | Front right (M1) | CCW | Servo 5 |
| 2 | Rear left (M2) | CCW | Servo 6 |
| 3 | Front left (M3) | CW | Servo 7 |
| 4 | Rear right (M4) | CW | Servo 8 |

Key ArduPilot parameters:
```
CAN_P1_DRIVER = 1
CAN_D1_PROTOCOL = 1
CAN_D1_UC_ESC_BM = 0x000F   ; bits 0–3 enable ESC nodes 1–4
Q_FRAME_CLASS = 1             ; Quad
Q_FRAME_TYPE = 1              ; X
```

Full parameter list in `docs/reference-guide/ardupilot-quadplane.md`.

#### Propellers

| Parameter | Value |
|---|---|
| Model | MAD FLUXER PRO 26×7.8 MATT — 2Pcs (SKU: 0103LX0013) |
| Diameter × pitch | 660.4 × 198.1 mm (26 in × 7.8 in) |
| Material | Carbon fiber + resin |
| Finish | Matte (non-reflective) |
| Weight | **69g per blade** — 4 installed props = **276g total** |
| Type | Fixed (not folding) |
| Hub bore | Ø10mm center shaft hole |
| Hub bolt pattern | Ø20mm bolt circle, 4× M3 screws |
| Accessories included | 8× M3×14mm mounting screws, 2× cover plates per pair |
| Optimum RPM | **1,900–4,000 RPM** |
| Single thrust limit | **18 kgf** — not the binding constraint (motor limits at 15 kgf) |
| Working temperature | -40°C to +65°C |
| Quantity ordered | 6 pairs (12 blades: 6 CW + 6 CCW); 4 installed, spares absorb early breakage |
| Unit price | $193.90/pair + $59.70 shipping |
| Reasoning (fixed not folding) | Folding props risk uncontrolled rotation in forward-flight airflow and uncertain locking behavior at speed |

**Hover operating point cross-check (prop test data vs motor datasheet):**

From the standalone prop test table:

| RPM | Thrust (gf) | Output Power (W) | Prop efficiency (gf/W) |
|---|---|---|---|
| 3,226 | 4,110 | 444 | 9.3 |
| 3,619 | 5,120 | 625 | 8.2 |
| 3,971 | 6,160 | 834 | 7.4 |
| 4,232 | 6,940 | 1,008 | 6.9 |

At hover (~3,780 RPM, 44.4V nominal): interpolated ~5,700 gf, prop efficiency ~7.8 gf/W. Motor datasheet showed 6,440 gf at 3,811 RPM at 48V — consistent after accounting for voltage difference. Numbers agree across both datasheets.

Hover operating point (~3,780 RPM) falls within the 1,900–4,000 RPM optimum range.

**Hub compatibility:** The prop hub requires a 10mm shaft bore. Confirm V8013 PRO output shaft diameter is 10mm before installation. Torque plate bolt pattern (Ø20mm, M3×4) must match the motor's prop adapter plate — verify with Alperen's motor mount design.

**Weight upgrade path:** Ultra-light carbon variants (~20g savings per prop) deferred to ~August 2026 after initial testing (80g total saving across 4 props).

---

### 4.6 HV Wiring Harness

#### Wire Sizing

Current ratings below assume silicone-insulated stranded copper in open-air UAV duty (not bundled in conduit).

| Segment | AWG | Est. Length | Est. Current | Connector |
|---|---|---|---|---|
| Battery main leads (positive + negative) | 4 AWG | 0.2–0.4 m | 320A peak, ~98–124A hover | AS150 or XT150 at battery; busbars at distribution end |
| HV bus to each ESC (×4) | **12 AWG** | up to 0.8 m (to boom root) | ~25–31A hover, 80A rated peak | XT90S or solder direct to ESC power leads |
| ESC power leads (factory, ESC to bus) | **12 AWG** | **800mm** (factory) | 80A rated | Bare tinned — add XT90S or solder to bus |
| ESC phase wires to motor (factory) | **14 AWG** | **150mm** (factory) | 80A rated (short duty only) | Bare tinned — 5mm bullet at motor side |
| BEC input from HV bus | 16 AWG | 0.2–0.3 m | 5A max | XT30 |

**Wire notes:**

The AMPX 80A ESC ships with 12AWG/800mm power leads and 14AWG/150mm phase leads. These are MAD's designed wire lengths for their thermal model at rated current.

- **Phase leads (14AWG/150mm): do not extend** without upsizing to 12AWG or heavier. The 150mm is short by design — heat dissipates through the motor housing. A longer 14AWG run at 80A will overheat.
- **Power leads (12AWG/800mm):** adequate for Spearhead's hover current (~25–31A/motor). At rated 80A the wires run warm but duty cycle is short (10–30s per motor zone limits). If the distribution bus is within 800mm of the ESC, terminate directly. If the boom run requires more length, splice with 10AWG.
- **Main battery leads (4AWG):** adequate for 200A continuous. At 320A peak (all 4 ESCs at rated max simultaneously, ~5s) thermal margin is acceptable for prototype duty. Monitor during first high-throttle runs. Upsize to 2AWG for sustained full-power testing.

Current rating reference (silicone insulated, open air):
- 4 AWG: ~200A continuous
- 10 AWG: ~55–65A continuous
- 12 AWG: ~40–50A continuous
- 14 AWG: ~32–40A continuous

#### CAN Bus Wiring

| Parameter | Value |
|---|---|
| Cable type | Twisted pair, 24 AWG, ≥ 10 twists/foot |
| Shielding | Single drain wire, grounded at Pixhawk end only |
| Connector | 4-pin JST-GH (CAN H, CAN L, 5V, GND) per DroneCode standard |
| Total run | ~2.5 m nose to tail — no signal integrity issue at 1 Mbps |
| Topology | Daisy-chain (node-to-node) with 120Ω terminators at each end of bus |
| Mid-bus nodes | No termination on intermediate nodes (ESCs, nose GPS) |

#### Servo Wiring (Tail Harness)

| Parameter | Value |
|---|---|
| Signal wire | 26 AWG, < 0.5 m per servo within tail bay |
| Power wire to tail | 18 AWG (7.8V servo rail, runs from BEC in nose bay to tail via harness) |
| Signal logic level | 5V (from Here4 breakout or tail PCB) |
| Power logic | 7.8V (separate supply, shared ground) |

#### Wing Panel Connector

The detachable wing panels each carry one flaperon servo. The connector at the wing root must break out signal, servo power, and ground.

| Signal | Notes |
|---|---|
| Servo signal (×1 per wing) | 5V logic PWM from nose FC |
| 7.8V servo power | From nose 7.8V BEC via fuselage harness |
| Ground | Shared with avionics ground |

Connector candidate: Molex Mini-Fit Jr 6-pin (known from Quiver, rated 9A per pin, latching). TBD in WP-E06.

---

### 4.7 Flight Controller

| Parameter | Value |
|---|---|
| Model | Pixhawk 6C (Holybro) |
| Firmware | ArduPlane 4.x with QuadPlane |
| IMU | Triple redundant (ICM-42688-P × 2, ICM-20649) |
| Barometer | Internal (ICP-20100) + Here4 external baro via CAN |
| Magnetometer | External via Here4 on CAN (internal mag disabled) |
| CAN ports | FMU-CAN1 (active), FMU-CAN2 (reserved) |
| UART ports | 6× (TELEM1, TELEM2, GPS1, GPS2, SERIAL4, SERIAL5) |
| PWM outputs | 16× (8× FMU direct, 8× via IO processor) |
| Power input | Dual: POWER1 (primary BEC/PM), POWER2 (backup BEC) |
| Vibration isolation | Internal gel-foam isolators in chassis |
| Location | Nose avionics bay, oriented with X-axis forward, mounted to isolator plate |
| USB-C | Ground access for configuration and parameter upload |

Power module: Holybro PM07 or PM06 V2 (12S input rated, 60A HV-side, provides 5.2V/6A to FC, voltage and current telemetry on ADC). Acts as POWER1 input and HV-side current sensor.

---

### 4.8 GPS and Navigation Sensors

#### Here4 GPS (×2)

| Parameter | Value |
|---|---|
| Model | CubePilot Here4 |
| GNSS constellations | GPS, GLONASS, Galileo, BeiDou |
| RTK capable | Yes |
| Interface | DroneCAN (primary), UART (secondary) |
| Internal sensors | Barometer, IMU, magnetometer |
| PWM outputs | Available on second cable via breakout board — used for tail CAN-to-PWM (see §4.10) |
| Supply voltage | 5V |
| CAN node IDs | Nose: 10, Tail: 11 |
| Separation | ~2 m (nose to tail) for GPS yaw accuracy (moving baseline or GPS-for-yaw mode) |

ArduPilot parameters:
```
GPS_TYPE = 9    ; DroneCAN
GPS_TYPE2 = 9
GPS_GNSS_MODE = 0     ; auto
GPS_GNSS_MODE2 = 0
GPS_AUTO_SWITCH = 1   ; auto-switch to best GPS
```

GPS yaw (moving baseline, dual-antenna heading) reduces reliance on magnetometer. Configure via `GPS_MB_UART_PORT` or CAN GPS heading — confirm supported with ArduPilot 4.x on DroneCAN Here4.

#### Airspeed Sensor

| Parameter | Value |
|---|---|
| Hardware | Pitot tube from storage (model TBD — Alperen to confirm) |
| Interface | I2C or analog ADC depending on pressure transducer |
| Location | Nose or boom, away from prop wash and wing downwash |
| ArduPilot | ARSPD_TYPE = 1 (analog) or 7 (I2C DLVR) depending on transducer |

If the stored pitot's transducer is unknown or incompatible, replacement option: Matek AP_ANALOG_AIRSPEED or Holybro airspeed sensor (I2C DLVR, ~$25).

#### Radar or LiDAR Altimeter

Selection (OQ-03):

| Model | Type | Range | Update rate | Interface | Weight | Cost |
|---|---|---|---|---|---|---|
| Ainstein US-D1 | Radar | 0.5–50 m | 12 Hz | UART | 56g | ~$250 |
| Benewake TF03-180 | LiDAR | 0.1–180 m | 100 Hz | UART | 76g | ~$250 |

US-D1: robust to dust, rain, and low-light. Better for operational use.
TF03-180: higher update rate and longer range. Better for precision hover.

Both integrate with ArduPilot via SERIAL_PROTOCOL = 9 (Lidar/rangefinder) or similar. Decision in WP-E05.

---

### 4.9 CAN Bus Architecture

CAN1 bus node map:

| Node ID | Device | Type | Notes |
|---|---|---|---|
| 1 | ESC1 — Motor 1 (FR) | AMPX 80A | DroneCAN ESCStatus + RawCommand |
| 2 | ESC2 — Motor 2 (RL) | AMPX 80A | |
| 3 | ESC3 — Motor 3 (FL) | AMPX 80A | |
| 4 | ESC4 — Motor 4 (RR) | AMPX 80A | |
| 10 | Here4 #1 (nose) | GPS + baro + compass | Primary GPS |
| 11 | Here4 #2 (tail) | GPS + baro + PWM out | Secondary GPS + tail CAN-to-PWM |

Bus termination: 120Ω at Pixhawk CAN port (built-in or terminator plug) and 120Ω at the tail Here4 (last node). No termination at intermediate nodes.

CAN2 is reserved. Future uses: generator controller (Phase 2), companion computer CAN bridge, additional sensors.

ESC node IDs are assigned via DroneCAN UI Tool (pre-flight on bench, not in the field). Set before boom installation.

---

### 4.10 Tail Electronics

**Problem:** Running PWM servo signal cables 2.5 m from the nose avionics bay to the tail surfaces through an environment with 4× 80A ESC switching noise causes signal degradation and potential interference.

**Solution:** CAN to the tail, PWM local to the tail. Here4 #2 serves as the tail CAN node and also sources PWM signals for V-tail servos via its second output cable.

#### Here4 CAN-to-PWM Path

The Here4 has two cables exiting the case:
1. CAN + voltage (currently connected on Quiver)
2. PWM/UART breakout (currently unused on Quiver)

Connecting both cables to the breakout board exposes PWM output pins. No case modification required.

**Validation test:** Thomas brings his Here4 to Ankara (~May 12). Zeynep's test procedure: connect Here4 to FC over CAN, connect Here4 second cable to breakout, confirm PWM signal outputs with scope or servo response. Post results to Spearhead forum. This test gates WP-E03.

#### Tail PCB (Erick — WP-E03)

A custom 2-layer PCB integrates all tail-side electrical functions:

| Function | Details |
|---|---|
| CAN input | 4-pin JST-GH, daisy-chained from CAN1 bus |
| 120Ω CAN termination | Jumper-selectable (end-of-bus position) |
| HV power input | HV bus spur (36–50V) via XT30 from tail boom harness. No regulated LV runs from nose to tail. |
| 7.8V regulation | HV → 7.8V switching regulator (e.g. LM5175 or equivalent module), ~3A capacity for 2 ruddervator servos |
| 5V regulation | 7.8V → 5V LDO or small switcher, ~500mA for Here4 |
| PWM breakout | 4× 3-pin servo headers (S / V+ / GND), signal at 5V from Here4 |
| Power connector (input) | XT30 — HV spur from rear boom ESC junction or dedicated tail harness branch |
| Servo connectors (output) | 4× standard 2.54mm 3-pin headers |
| Board size target | ~60 × 40 mm (slightly larger to accommodate switching regulator) |
| Layers | 2 |

Voltage hierarchy: HV in → 7.8V (ruddervator servos) → 5V (Here4). All three rails share a common ground. PWM signal from Here4 at 5V logic is compatible with servos powered at 7.8V with shared ground (confirmed by Erick).

**Fallback Option A (if Here4 test fails): Matek CAN-L4-RC**

Matek CAN-to-PWM or SBUS converter board, ~$35. Connect CAN bus to this board and run PWM outputs to servos. No custom PCB required. Slightly less integrated but eliminates dependency on Here4 PWM cable.

**Fallback Option B (first hover test only):** Point-to-point wiring from Here4 breakout to servos, with BECs hanging loose in the tail bay. Servos will be slightly underpowered (connected 5V rail via Here4, not optimal 7.8V rail) but functional for basic hover validation. Replace with tail PCB before any aggressive maneuvering.

Decision checkpoint: after Thomas's Here4 test (~May 12). If test passes, proceed with tail PCB (WP-E03). If not, implement Matek board.

---

### 4.11 Servo System

#### Ruddervator Servos (×2)

V-tail airfoil: NACA 0015, 280mm chord, control surface 30–35% chord (~84–98 mm deep surface). Each surface is one ruddervator providing mixed pitch + yaw authority (ArduPilot SERVO_FUNCTION 79/80).

| Parameter | Requirement |
|---|---|
| Type | Digital HV (7.4–8.4V operating range) |
| Torque target | ≥ 8 kg·cm at 7.4V |
| Speed | ≤ 0.12 s / 60° at 7.4V |
| Weight | ≤ 30g each (tail weight is CG-critical) |
| Connector | Standard 3-pin JR/Futaba (2.54mm pitch) |
| Control | Via tail PCB PWM output (see §4.10) |

Torque target is preliminary — Zeynep's control surface sizing study output (WP-E06 input) will set the final requirement. Verify hinge moment calculation before locking servo selection.

#### Flaperon Servos (×2)

Main wing: Clark Z, 20–25% chord flaperons.

| Parameter | Requirement |
|---|---|
| Type | Digital HV |
| Torque target | ≥ 6 kg·cm at 7.4V |
| Speed | ≤ 0.15 s / 60° at 7.4V |
| Weight | ≤ 30g each |
| Connector | Standard 3-pin + wing panel connector (see §4.6) |
| Control | Via Pixhawk IO PWM outputs 1 and 2 |

Servo candidates (both types): Savöx SH-0255MG (~20g, 9 kg·cm at 7.4V), KST DS135MG (~27g, 13 kg·cm at 7.4V), MKS HV6130 (~20g, 10 kg·cm at 7.4V). Selection in WP-E06.

---

### 4.12 Telemetry and Communications

#### RC Link

| Parameter | Value |
|---|---|
| Range requirement (Phase 1) | VLOS, ≤ 5 km |
| Protocol | CRSF (preferred, bidirectional) or SBUS |
| FC interface | UART5 (SERIAL5) |
| Receiver location | Nose bay or fuselage side, antenna external |
| Candidates | ExpressLRS (900 MHz, CRSF, 1W output option) or TBS Crossfire |

ExpressLRS preferred: native CRSF bidirectional allows MAVLink passthrough to GCS without a separate telemetry radio at short range. CRSF also gives better latency than SBUS. Confirm transmitter compatibility with team.

#### Telemetry Radio

| Parameter | Value |
|---|---|
| Range (Phase 1) | 30–80 km LOS with whip antennas |
| Protocol | MAVLink 2.0 |
| FC interface | UART1 (TELEM1) |
| Data | Attitude, GPS, battery SoC, ESC temps, link quality |
| Candidates | RFD900x (1W, 900 MHz, ~$200/pair), SiK 915 MHz (100mW, ~$80/pair) |

For Phase 1 VLOS testing, SiK 915 MHz is adequate and lower cost. Upgrade to RFD900x before Phase 3 extended-range runs (OQ-04 / WP-E08).

Phase 3 requirement (> 200 km) likely requires satellite. See §6.1.

---

### 4.13 Health Monitoring

#### Battery Monitoring

Smart battery over DroneCAN — no external power module or ADC wiring required.

| Parameter | Method |
|---|---|
| Pack voltage | Smart battery BatteryInfo → DroneCAN → FC |
| Pack current | Smart battery BatteryInfo → DroneCAN → FC |
| Cell voltages | Smart battery BatteryInfo → DroneCAN → FC (per-cell monitoring, not available with a standard PM07) |
| SoC | Reported by smart battery BMS (fuel gauge IC) — more accurate than coulomb counting |
| Temperature | Smart battery BMS sensor → DroneCAN → FC |
| Low voltage warning | BATT_LOW_VOLT = 42V (3.5V/cell), GCS alert + buzzer |
| Critical voltage | BATT_CRT_VOLT = 39.6V (3.3V/cell), RTL trigger |
| Low capacity | BATT_LOW_MAH = 2000 mAh (reserve threshold) |

ArduPilot parameter: `BATT_MONITOR = 8` (UAVCAN/DroneCAN). The smart battery node ID must be unique on the CAN bus. Confirm the selected unit broadcasts standard UAVCAN BatteryInfo or DroneCAN BatteryInfoAux messages — not all smart batteries use the standard message type.

#### ESC Telemetry and Thermal

AMPX 80A V2 DroneCAN ESCs broadcast ESCStatus messages in real time: voltage, current, temperature, and operation status per ESC. ArduPilot logs these automatically and streams over MAVLink telemetry. No additional wiring required.

**ESC thermal thresholds** (source: AMPX 80A ESC datasheet — these are ESC junction temperature limits, independent of motor winding temperature):

| ESC junction temp | Action |
|---|---|
| > 125°C | ESC reduces output to 50% max power; fault signal broadcast over DroneCAN |
| > 140°C | Full motor shutdown — will not restore until throttle zeroed AND ESC cools below 80°C |

ESC ambient operating limit is 65°C. In warm weather (Ankara summer) with ESCs mounted externally on booms, ensure the heatsink fins have unobstructed airflow during hover. Log ESC temperature via DroneCAN telemetry and configure GCS alerts well below the 125°C threshold (suggest 90°C warning).

**CAN watchdog — critical:** The AMPX 80A cuts motor power if CAN commands are absent for more than **2 seconds**. A CAN bus fault mid-hover causes all four motors to cut simultaneously after 2s. Verify `CAN_D1_UC_ESC_BM` is correct and confirm CAN bus health before first hover.

#### Motor Thermal

The V8013 PRO motor has no onboard temperature sensor and no direct telemetry output for winding temperature. The operating zone limits below come from the **motor datasheet footnote** (not the ESC) and describe the motor's sustainable current based on winding copper losses and prop-wash cooling at 25°C ambient:

| Motor current | Zone | Duration |
|---|---|---|
| < 24A/motor | Sustainable continuous | Indefinite |
| 24–78A/motor | Short-term | ~10–30s |
| > 78A/motor | Non-working | Avoid |

**These limits are independent of the ESC.** The ESC will run at 80A continuous without complaint; it is the motor windings that constrain sustained current.

As shown in §4.1, the motor runs in the continuous zone (< 24A) for most of the useful flight window and only approaches short-term territory near the low-battery warning voltage. The zones are more of a concern at altitude, where throttle must increase to compensate for lower air density.

**Monitoring approach (no direct motor temp sensor):**

| Proxy | How | Limitation |
|---|---|---|
| Current per motor | ESC DroneCAN telemetry (ESCStatus.current) | Best real-time indicator of which zone you're in |
| Motor case temperature | Hand check after landing (bench tests only) | Rudimentary but effective for initial characterization |
| RPM | ESC DroneCAN telemetry | Confirms expected hover RPM vs. throttle; anomaly = mechanical issue |

For Phase 2 characterization: consider adding NTC thermistors to the motor stator cans, routed along the boom wiring to the FC ADC, to establish actual winding temperature margins before committing to operational flight profiles.

#### Pre-flight Diagnostics

ArduPilot built-in pre-arm checks cover: GPS lock and HDOP threshold, compass calibration, RC calibration and failsafe, ESC arming, battery voltage above BATT_ARM_VOLT, CAN node health, and altimeter sanity check. No additional software needed for Phase 1.

#### Phase 2 Additions (stub)

Engine temp: type-K thermocouple at cylinder head and exhaust, via MAX6675 or AD8495 amplifier to FC ADC. Fuel level: ultrasonic or float sensor to ADC or UART.

---

### 4.14 Heading LEDs

| Parameter | Value |
|---|---|
| Type | WS2812B (NeoPixel) addressable RGB |
| Control | ArduPilot NTF (Notify) subsystem |
| FC interface | GPIO pin configured as NeoPixel output |
| ArduPilot parameter | NTF_LED_TYPES = 1 (NeoPixel) |
| Color scheme | Front/port: red, front/starboard: green, rear: white/flashing red |
| Location | Wing tips (port/starboard), nose, tail |
| Power | 5V avionics rail; WS2812B draws ~60 mA per LED at full white |
| Estimated LEDs | 4–6 per wing tip, 2 at nose, 2 at tail |
| Total power budget | ~16 LEDs × 20 mA average = ~320 mA at 5V |

ArduPilot NTF patterns: arming flash, mode color (VTOL = yellow, FW = blue, RTL = red), GPS quality indicator.

---

### 4.15 Component Locations

Detailed layout TBD in WP-E07, subject to Alperen's fuselage truss dimensions. Preliminary placement:

```
Nose bay (forward of CG):
  - Pixhawk 6C: center, vibration-isolated, X-axis forward
  - Here4 #1 GPS: top of nose, unobstructed sky view
  - RC receiver: side-mounted, antenna external
  - Telemetry radio: side-mounted, antenna external
  - 5V BEC (×2) + 7.8V BEC: floor of nose bay
  - Power module (PM07): inline with battery lead entry point

CG bay (center, under wing):
  - 12S battery: centered fore-aft and laterally on CG
  - HV power bus / distribution PCB: adjacent to battery
  - HV kill switch: accessible from exterior hatch

Boom roots (×4):
  - ESCs: as close to motors as boom structure allows (reduces phase wire length)
  - Phase wire exits sealed where they enter boom
  - CAN bus daisy-chained along boom exterior or inside boom if square tube allows

Tail bay:
  - Here4 #2 GPS: top of tail, unobstructed sky view
  - Tail PCB: mounted to bulkhead, connectors accessible

Wings (detachable panels):
  - Flaperon servo: inline with control horn, servo output axis on hinge line
  - Wing panel connector: at wing root, mates to fuselage connector on insertion
```

---

### 4.16 Phase 1 BOM and Cost Estimate

**Already ordered** (status: in transit, MAD Motor Poland, ~May 1, 2026):

| Item | Qty | Est. Unit Cost | Est. Total |
|---|---|---|---|
| MAD V8013 PRO IPE 150KV motor | 6 | ~$185 | ~$1,110 |
| MAD AMPX 80A DroneCAN ESC (V2) | 6 | $60 | $360 |
| MAD FLUXER PRO 26×7.8 MATT (per pair) | 6 pairs | $193.90 | $1,163 + $60 shipping |
| **Subtotal (ordered)** | | | **~$2,633** |

**To procure:**

| Item | Qty | Est. Unit Cost | Est. Total |
|---|---|---|---|
| Pixhawk 6C (Holybro, no PM07 needed) | 1 | ~$360 | ~$360 |
| 12S DroneCAN smart battery ≥ 22Ah | 1 | TBD | TBD |
| GPS module — nose (TBD: Here4 or alternative) | 1 | TBD | TBD |
| Here4 GPS — tail (CAN-to-PWM + GPS2) | 1 | $300 | $300 |
| 5V/10A switching BEC (×2 for redundancy) | 2 | $25 | $50 |
| 7.8V/10A switching BEC (flaperon servos, CG bay) | 1 | $30 | $30 |
| HV kill switch + contactor (or PP75) | 1 | $50 | $50 |
| LV kill switch | 1 | $15 | $15 |
| AS150 connector pairs | 2 | $12 | $24 |
| XT90S connector pairs (ESC boom junctions + spares) | 8 | $6 | $48 |
| 5mm bullet connectors (×12 sets) | 12 | $2 | $24 |
| 4 AWG silicone wire (2 m) | 1 lot | $15 | $15 |
| 8 AWG silicone wire (10 m) | 1 lot | $20 | $20 |
| 12 AWG silicone wire (extension leads) | 1 lot | $12 | $12 |
| LV wiring assortment (18–26 AWG) | 1 lot | $30 | $30 |
| CAN bus twisted-pair wire (5 m) | 1 lot | $10 | $10 |
| JST-GH crimp kit (CAN, UART, assorted) | 1 | $25 | $25 |
| Digital HV ruddervator servos | 2 | $40 | $80 |
| Digital HV flaperon servos | 2 | $35 | $70 |
| Wing panel servo connectors (Molex Mini-Fit Jr) | 2 sets | $10 | $20 |
| RC receiver (TBD — pending RC system decision) | 1 | TBD | TBD |
| Telemetry radio (TBD — pending SIYI HM30 decision) | 1 | TBD | TBD |
| Radar or LiDAR altimeter (TBD) | 1 | ~$280 | ~$280 |
| Pitot tube transducer (if needed beyond existing) | 1 | $30 | $30 |
| WS2812B LED strip + wiring | 1 lot | $20 | $20 |
| Tail PCB fabrication (JLCPCB, 5 boards) | 5 | $2 | $10 |
| Tail PCB component BOM | 1 lot | $30 | $30 |
| Matek CAN2PWM (fallback, keep in reserve) | 1 | $35 | $35 |
| Heat-shrink, zip ties, mounting hardware | 1 lot | $40 | $40 |
| **Subtotal (to procure, est.)** | | | **~$1,578 + TBD items** |
| **Phase 1 Grand Total (est., excl. labor)** | | | **~$4,211 + TBD items** |

Estimates only. Verify against actual quotes before procurement. Some items may be on hand from Quiver or other projects.

---

## 5. Phase 2: Hybrid Integration

Target: July–October 2026. Objectives: IC engine integration, transition flights, generator architecture decision and initial implementation, vibration characterization and notch filter tuning.

> **⚠️ STATUS: This section is preliminary. Significant research is still required before any of the Phase 2 content below should be treated as a design baseline.**
>
> The IC engine, fuel system, starter, generator, ignition, RPM sensor, throttle servo, and vibration isolation subsections in §5 are working notes based on general 2-stroke RC engine practice and a candidate engine class (35–55cc Stinger / DLE). The actual engine has not been selected or purchased, no bench testing has occurred, and most parameters (fuel consumption rate, electrical load profile, vibration spectrum, alternator output, starter inrush current) will only be known after empirical measurement on the chosen unit. Treat all numbers, parts, and architectures in §5 as candidate values to be validated, not commitments. Specific items requiring further research before commitment:
>
> - Engine selection (Stinger vs DLE, 35cc vs 55cc) — buy both and test empirically per current plan
> - Starter motor selection and current draw (Pilot RC vendor response pending)
> - VESC dual-role starter-generator viability (no off-the-shelf precedent at this scale)
> - Generator output topology (alternator vs starter-as-generator) and rectifier specifications
> - Ignition battery: shared 12S vs dedicated pack
> - CHT / EGT sensor selection and ArduPilot integration
> - Real vibration spectrum and final notch filter tuning (cannot be predicted, must be measured)
> - Throttle servo torque and stroke for the actual carburetor lever geometry
>
> Phase 1 (§4) is the working design baseline. Phase 2 will be re-scoped after Phase 1 hover testing and once the engine candidates are in hand.

---

### 5.1 Fuel System and Mix

#### 2-Stroke Fuel Mix

DLE and Stinger engines in the 35–55cc class are air-cooled single-cylinder 2-stroke gasoline engines. They require a premixed fuel-oil mixture — no separate oil reservoir.

| Stage | Ratio (fuel:oil) | Notes |
|---|---|---|
| Break-in (first 2–3 tanks) | 25:1 | Richer oil mix lubricates and seats rings; run at varied throttle, avoid sustained WOT |
| Running | 40:1 to 50:1 | Manufacturer-specific; DLE 35cc manual specifies 40:1 typical |

**Oil type:** semi-synthetic or full-synthetic 2-stroke oil designed for air-cooled engines. Do not use 4-stroke oil, marine 2-stroke (designed for water-cooled), or lawnmower 2-stroke oil. Brands: Motul 800 2T, Castrol Power 1 Racing 2T, or equivalent.

**Fuel:** 92+ octane unleaded gasoline, ethanol-free preferred. Ethanol is hygroscopic (absorbs atmospheric moisture), corrodes aluminum carburetors and fuel lines over time, and degrades stored fuel faster. Use ethanol-free pump gas (E0) or aviation mogas where available.

**Altitude correction:** 2-stroke carburetors are calibrated for sea level. At 1,200m operational altitude (air density ~88% of sea level), the mixture runs slightly rich at the same needle position. Most engines have a main jet needle clip with 3–5 positions for mixture adjustment. Research required: determine the DLE 35cc / Stinger jet needle adjustment procedure and characterize the correction needed at operational altitude. An engine running rich at altitude produces less power, runs hot in a different way, and fouls plugs.

**Pre-flight mixing protocol:** Mix fuel fresh before each test campaign. Fuel-oil premix has a shelf life of approximately 30 days before oil separation begins. Do not use old premix.

#### Fuel Storage and Delivery

| Item | Notes |
|---|---|
| Tank material | Aluminum or HDPE (compatible with gasoline). Oratex and standard foam are NOT compatible with gasoline — design fuel bay accordingly. |
| Tank volume | Target ~6L (per Alperen's fuselage layout). At ~35cc/min consumption estimate for DLE 35 at cruise, 6L gives ~170 min endurance (needs empirical verification). |
| Fuel lines | Viton or Tygon fuel tubing (gasoline-compatible). Standard silicone tubing is NOT gasoline-compatible and will degrade. |
| Fuel filter | Inline filter between tank and carb. Replace every 5–10 flights. |
| Clunk pickup | Weighted clunk at tank inlet ensures fuel feed at any aircraft attitude. |
| Vent | One-way vented cap to prevent vapor lock. |

---

### 5.2 Throttle Control

The throttle arm on the engine carburetor is driven by a servo, controlled via ArduPilot.

| Parameter | Value |
|---|---|
| Servo type | Digital, metal gear, standard torque (≥ 3 kg-cm) |
| ArduPilot function | `SERVO9_FUNCTION = 70` (throttle) — maps ICE throttle to servo output 9 |
| Servo travel | Must be calibrated to match carburetor idle stop (minimum) and wide-open throttle (maximum). Set `ICE_IDLE_RPM` and test with tachometer. |
| Idle RPM | Engine-specific; DLE 35cc target ~1,400–1,600 RPM at idle. Low enough to not produce significant thrust, high enough to stay running. |
| Throttle during VTOL hover | ArduPilot keeps IC engine at idle during all VTOL modes (`ICE_IDLE_RPM` via `ICE_ENABLE`). Engine is not shut down during hover — it idles at cruise throttle stop. |
| Throttle during transition | ArduPilot advances IC throttle as VTOL motors ramp down across `Q_TRANSITION_MS`. |
| Emergency cutoff | Ignition kill (see §5.3). Throttle servo drives to idle position as secondary measure. |

ArduPilot IC engine parameters (research required — no direct experience with ICE ArduPilot integration):
```
ICE_ENABLE = 1
ICE_START_CHAN = <RC channel for starter relay>
ICE_IDLE_PCT = <idle throttle %, calibrate empirically>
ICE_IDLE_RPM = 1500   ; target idle RPM, uses RPM sensor feedback
ICE_IDLE_DB = 50      ; RPM deadband ±50 RPM
ICE_RPM_CHAN = 1      ; RPM sensor instance
```

Confirm ArduPilot ICE library supports DroneCAN ESC coexistence — the ICE library was designed for planes with a single IC engine and may need validation in a QuadPlane configuration.

---

### 5.3 Choke and Cold-Start Sequence

Most DLE and Stinger engines use a manual choke plate (butterfly valve) in the carburetor air intake. There is no electronic choke — it must be mechanically actuated.

#### Choke Actuation Options

| Option | Description | Pros | Cons |
|---|---|---|---|
| Manual lever | Ground crew or pilot operates choke arm directly before flight | Zero weight/complexity | Requires physical access to engine; impractical after integration in fuselage |
| RC-servo actuated | Micro servo on choke arm, assigned to a transmitter switch | Remote operation, no hardware complexity | Adds 1 servo; needs physical routing inside engine bay |
| Pull cable | Bowden cable from cockpit hatch to choke arm | Simple, no electronics | Cable routing through fuselage, potential binding |

**Recommendation for Phase 2 prototype:** servo-actuated choke on a dedicated RC channel (transmitter 2-position switch). Assign as a passthrough output, not an ArduPilot-controlled function. Pilot operates it manually during the start sequence.

#### Cold-Start Sequence (proposed, to be refined with engine testing)

1. Transmitter on, aircraft in QLOITER or QSTABILIZE, VTOL motors armed but at idle
2. Confirm throttle servo at idle position (minimum)
3. Close choke (RC switch → choke servo)
4. Engage electric starter (RC switch → starter relay) for 2–5s
5. Engine fires → immediately open choke partway
6. Run at fast idle (~2,000 RPM) for 30–60s to warm up
7. Open choke fully — engine now at normal idle
8. Confirm RPM sensor reading matches expected idle RPM
9. Proceed with VTOL takeoff as normal

**Hot restart (engine already warm):** skip choke, engage starter briefly.

Note: If the engine does not start within 3 attempts, investigate flooding (excess fuel in cylinder). Recovery: open choke fully, wide-open throttle, crank briefly to clear — a procedure that should only be done on the ground with props removed.

---

### 5.4 Ignition System

2-stroke engines at this scale use a CDI (Capacitor Discharge Ignition) module mounted on the engine. The CDI generates high-voltage spark from a low-voltage input and a trigger from a hall effect sensor or pickup coil on the flywheel.

| Item | Notes |
|---|---|
| Power supply | CDI typically requires 5–6V regulated DC, < 1A. Use a dedicated low-noise 5V BEC (separate from the avionics BEC to avoid back-EMF transients on the avionics rail). |
| Kill wire | Most CDI modules have a "kill" wire: ground this wire to stop ignition immediately. Connect to an FC relay output (`SERVO10_FUNCTION = 28` relay) and an accessible hardware switch. Both the pilot RC switch and the hardware switch must independently kill the engine. |
| Spark plug | Replace per engine manual (typically every 10–20 flight hours). Inspect for fouling during initial testing. |
| RF interference | CDI spark generates broadband RF noise. Shield ignition wires with braided metal sleeve, ground the shield at one end. Route CDI wires away from GPS antennas and CAN bus cables. Confirm GPS fix quality doesn't degrade with engine running during ground test. |

---

### 5.5 Electric Starter

#### Pilot RC Auto-Starter (primary candidate)

The Pilot RC DLE35RA electric auto-start unit (420g) is a gear-driven brushless motor that engages the engine crankshaft for starting. Current status: vendor email sent, no response as of May 2026. May be discontinued or small-batch.

**Starter power circuit:**

```
[2S LiPo, 7.4V, ~1,500–2,000mAh]
        │
  [High-current relay, 30–50A]  ←── FC relay output (SERVO10 or AUX)
        │                              or dedicated RC transmitter switch
  [Starter motor]
        │
  [Engine crankshaft]
```

| Parameter | Value |
|---|---|
| Supply voltage | 7.4V (2S LiPo) or 11.1V (3S) — confirm with Pilot RC spec |
| Peak current | ~20–40A for < 3s (estimated; measure at bench) |
| Relay | 30–50A relay (Axiom, Tyco, or automotive relay). Switching inductive load — add flyback diode across relay coil. |
| FC control | ArduPilot `ICE_START_CHAN` assigns an RC channel. When triggered, the FC engages the relay, waits for RPM sensor to confirm start, then disengages. |
| Battery isolation | Dedicated 2S or 3S LiPo prevents starter inrush from sagging the main 12S VTOL battery. Also electrically isolates a potential starter motor failure from the avionics rail. |
| Safety interlock | Starter must only be able to engage when engine throttle servo is at idle and aircraft is on the ground (or explicitly commanded in-flight for in-flight restart). Use ArduPilot `ICE_STARTTIME` parameter to limit maximum starter engagement duration. |

#### Fallback: Manual Pull-Start

If the Pilot RC unit is unavailable, the engine is pull-started on the ground before arming. The aircraft must not be armed (VTOL motors must not be able to spin) during pull-starting. This is standard practice for Phase 2 initial testing regardless of whether an electric starter is installed.

---

### 5.6 RPM Sensor

The RPM sensor is a **Phase 2 day-one requirement**, not optional. ArduPilot uses it for:
1. `ICE_IDLE_RPM` closed-loop idle control
2. `INS_HNTCH_MODE = 3` — harmonic notch filter uses RPM to dynamically track engine vibration frequency
3. Engine health monitoring and telemetry
4. Fuel consumption estimation (if correlated with throttle position)

#### Sensor Options

| Option | Interface | How it works | Pros | Cons |
|---|---|---|---|---|
| CDI tach wire | FC GPIO, interrupt | Most CDI modules output a 1-pulse/rev signal on a tach wire | No added hardware, directly from engine | Confirm DLE/Stinger CDI has tach output; may require pull-up resistor |
| Hall effect sensor | FC GPIO, interrupt | Magnet epoxied to flywheel; hall sensor mounted on crankcase | Very reliable, clean signal | Requires access to flywheel for magnet installation |
| Optical sensor | FC GPIO, interrupt | Reflective tape on prop shaft; IR emitter/detector pair | Easy to install | Sensitive to oil/dirt fouling |

**Recommendation:** try CDI tach wire first (zero added hardware if it exists). If unavailable, use hall effect on the flywheel.

ArduPilot parameters:
```
RPM1_TYPE = 1         ; GPIO interrupt-based
RPM1_PIN = <GPIO>     ; pin number on Pixhawk 6C
RPM1_SCALING = 1.0    ; 1 pulse per revolution
RPM1_MAX = 10000      ; max expected RPM
RPM1_MIN = 500        ; min RPM to register as running
```

For harmonic notch filter:
```
INS_HNTCH_ENABLE = 1
INS_HNTCH_MODE = 3    ; tracks RPM sensor
INS_HNTCH_FREQ = 100  ; starting estimate (confirm with INS_VIBE log data)
INS_HNTCH_BW = 40
INS_HNTCH_ATT = 40
```

The notch filter must be tuned from real INS_VIBE log data taken during the first ground engine runs. Estimate: DLE 35cc at 6,000–8,000 RPM cruise → 100–133 Hz fundamental. Harmonics at 200–266 Hz may also require a second notch.

---

### 5.7 Engine Health Monitoring

| Sensor | Interface | ArduPilot integration |
|---|---|---|
| Cylinder head temperature (CHT) | K-type thermocouple → MAX6675 or MAX31855 SPI amplifier → FC SPI | Log via scripting or companion; no native ArduPilot CHT support — use Lua script to read SPI and forward via MAVLink |
| Exhaust gas temperature (EGT) | K-type thermocouple (higher range needed: 0–900°C) → MAX31855 | Same as CHT |
| Fuel level | Ultrasonic level sensor or float potentiometer → FC ADC | `BATT2_MONITOR = 7` (fuel flow, if using flow meter) or analog level to ADC for display |
| RPM | See §5.6 | Native ArduPilot RPM library |

**CHT limits (2-stroke, air-cooled, approximate):**
- Normal cruise: 150–200°C
- High: 200–230°C (acceptable for short time)
- Critical: > 250°C — engine damage risk; reduce throttle or enter VTOL mode

EGT is a better real-time mixture indicator than CHT. Rich mixture: lower EGT. Lean mixture: higher EGT and higher CHT. Target EGT range for 2-stroke at cruise: 550–680°C (engine-specific; characterize during break-in).

---

### 5.8 Generator Architecture and VESC Starter-Generator

#### Power Budget

| Load | Cruise power draw |
|---|---|
| Avionics + sensors | ~30–50W |
| Servos (cruise trim, minor inputs) | ~5–15W |
| VTOL battery top-up (trickle) | 50–150W (goal) |
| **Total generator target** | **100–200W** |

#### Architecture Decision: VESC Dual-Role Starter-Generator (Option A)

This is the primary architecture to develop. The same brushless motor and VESC controller handles both the electric start function and in-flight power generation, sharing shaft access to the IC engine output.

**Concept:**

```
[IC Engine crankshaft / prop shaft]
         │
   [Overrunning clutch or direct coupling]
         │
   [Brushless motor — starter/generator]
         │
   [VESC controller]
    ├── Starting mode: battery → VESC drives motor → spins engine
    └── Generator mode: engine drives motor → VESC rectifies/regulates → charges battery bus
```

**Motor and VESC selection:**

At engine cruise RPM of ~6,000–8,000 RPM, the generator motor must produce useful output voltage. For VESC to regulate a target bus voltage (say 50V for 12S charging), the motor KV should be selected so back-EMF at cruise RPM is in the right range.

At 7,000 RPM and target ~50V output (no-load):
- KV = RPM / V = 7,000 / 50 = **140 KV** (similar to the VTOL motors)
- At 7,000 RPM with 140 KV motor: VESC manages field weakening and duty cycle to regulate output

Alternatively, a lower-KV motor allows operation with lower DC bus voltage:
- 200–400 KV motor at 7,000 RPM: back-EMF ~17–35V; use a boost converter between VESC output and 12S battery

**VESC controller options:**

| Model | Cont. current | Voltage | Weight | CAN | Notes |
|---|---|---|---|---|---|
| VESC 6 MkVI | 50A cont | 60V max | ~100g | Yes | Well-supported; open source firmware; widely used in research |
| Flipsky VESC 6.7 | 80A cont | 60V max | ~80g | Yes | Lower cost VESC 6 derivative |
| VESC 75/300 | 300A cont | 90V max | ~300g | Yes | Overkill for this application |

For 100–200W generator: VESC 6 at ~50V/3–4A — well within 50A continuous rating.

**Shaft coupling and disengagement:**

The biggest unresolved mechanical question. Three options:

| Option | Description | Issues |
|---|---|---|
| Direct coupling, always engaged | Motor shaft permanently connected to engine output. VESC switches between motor and generator mode. | During VTOL hover, engine at idle (~1,500 RPM) drives the motor passively. VESC must be in freewheeling/passive mode to avoid drag. At engine shutdown, motor coasts freely. |
| Overrunning (sprag) clutch | One-way bearing. Motor drives engine for starting (motor faster than engine). Engine drives motor for generation (engine faster). | Sprag clutch direction must be correct. If motor spins faster than engine during starting, clutch engages. During cruise, engine spins faster, clutch engages other way. Need to verify this orientation works or use bidirectional clutch. |
| Electromagnetic clutch | Electric clutch that engages on command | Clean control; adds ~150–300g and complexity |

**Recommendation:** start with direct coupling and VESC freewheeling mode for Phase 2. The drag penalty at idle RPM is small (motor back-EMF is low at 1,500 RPM). If idle drag is unacceptable, add electromagnetic clutch in Phase 3.

**Battery charging from generator:**

The generator output does not use a conventional CC/CV LiPo charger. Instead:

```
VESC output (regulated DC, ~50V)
        │
  [Schottky diode, 5–10A, 60V]  ← prevents battery backfeed into VESC
        │
  [12S battery bus]
```

With the VESC regulating its output to slightly above the current battery bus voltage (e.g., bus at 44V → VESC regulates to 45V → ~1–3A flows into battery), the battery's own internal resistance limits the charge current naturally. This is a **parallel-source top-up**, not a dedicated charger. It is safe for in-flight use because:
- Current is self-limiting by voltage differential
- Battery is never charged above its full-charge voltage during flight (discharge is happening simultaneously)
- Schottky diode prevents reverse current if VESC output drops

VESC parameters for generator mode: set output voltage target via VESC Tool → `App Configuration → ADC` or via CAN commands from companion computer. Target: ~50V (slightly above nominal 44.4V hover bus voltage).

**Telemetry:** VESC supports CAN bus telemetry (voltage, current, RPM, temperature, power). Route to CAN2 on Pixhawk or to companion computer UART.

#### Alternative: Dedicated Belt-Driven Alternator (Option B)

A small permanent-magnet alternator (PMA) driven by a toothed belt from the engine prop shaft. Research needed:

- Candidates: small UAV PMAs from Genasun, Elroy Air, or custom-wound stators
- At 100–200W and 6,000–8,000 RPM, a 3-phase PMA produces AC → 3-phase bridge rectifier → voltage regulator → DC bus
- Weight target: < 300g for PMA + rectifier + regulator
- No engagement mechanism needed — belt slips or is tensioned by design
- Simpler than VESC but less flexible (fixed voltage-RPM curve unless regulator compensates)

Decision: develop Option A (VESC) in parallel with researching Option B. Select after bench testing both by early Phase 2.

---

### 5.9 Vibration Isolation and Notch Filter

**Mechanical isolation:**

- IC engine mounted on a soft sub-frame isolated from the main fuselage truss using silicone anti-vibration mounts (e.g., Lord Micro Mounts or M3 rubber standoffs)
- Do not hard-mount the engine to the carbon tube truss — transmitted vibration at 100–133 Hz will saturate the IMU
- Fuel tank should also be soft-mounted or suspended (liquid sloshing dampens vibration naturally)

**Electronic isolation:**

- Pixhawk 6C has triple-IMU redundancy with internal gel isolators. Adequate as a baseline.
- If INS_VIBE logs show clipping (values > 30 m/s²) during ground engine runs, add external foam pad isolation below the FC mount.
- Ignition CDI and spark plug wire: wrap in braided copper shielding, grounded at one end only. RF from ignition is broadband and will affect GPS if unshielded.
- Separate 5V BEC for CDI from avionics BEC. CDI draws pulsed current at spark frequency; this creates noise on shared supplies.

**Notch filter setup procedure:**

1. Engine running on bench (props OFF), aircraft secured
2. Log `INS_VIBE` and `INS_RAW_ACC1/2/3` at 50–100Hz for 60s
3. FFT of accelerometer data to identify vibration peaks
4. Set `INS_HNTCH_FREQ` to the primary peak frequency
5. Enable `INS_HNTCH_MODE = 3` (RPM sensor tracking)
6. Set `INS_HNTCH_HMNCS = 3` (fundamental + first harmonic)
7. Recheck `INS_VIBE` with notch enabled — peaks should attenuate by > 20 dB
8. Log comparison: before/after notch, at multiple throttle positions (idle, 50%, WOT)

---

### 5.10 Phase 2 Preliminary BOM

| Item | Est. Cost | Notes |
|---|---|---|
| IC engine — DLE 35cc (or Stinger 35cc) | ~$200 | Buy two sizes: 35cc and 55cc for comparison testing |
| IC engine — DLE/Stinger 55cc | ~$280 | |
| Carburetor jet needle set (altitude adjustment) | ~$20 | For altitude mixture correction research |
| Gasoline (E0, 92+) | ~$10/gal | Ongoing cost |
| 2-stroke synthetic oil (Motul 800 or equivalent) | ~$25/L | One liter covers ~40L of fuel at 40:1 |
| Ignition CDI + kill wire wiring | ~$50 | Usually included with engine |
| Throttle servo (metal gear digital, small) | ~$20 | |
| Choke servo (micro, < 10g) | ~$15 | |
| Choke servo mount | ~$10 | 3D printed |
| High-current relay (30–50A) + flyback diode | ~$15 | For starter control |
| Electric starter (Pilot RC DLE35RA) | ~$250 | If vendor responds; otherwise pull-start for Phase 2 |
| Starter 2S LiPo (1500–2000mAh) | ~$20 | Isolated starter supply |
| K-type thermocouple (×2: CHT + EGT) | ~$10 | |
| MAX31855 SPI thermocouple amplifier (×2) | ~$15 | |
| RPM hall effect sensor or optical sensor | ~$10 | If CDI tach output unavailable |
| Inline fuel filter | ~$5 | |
| Viton fuel tubing (1m) | ~$10 | Gasoline-compatible |
| Fuel tank (6L, aluminum or HDPE) | ~$40 | |
| Engine vibration isolator mounts (Lord or silicone) | ~$30 | |
| Separate 5V BEC for CDI | ~$15 | Noise isolation from avionics |
| VESC 6 MkVI or Flipsky 6.7 | ~$120 | Starter-generator controller |
| Brushless motor for starter-generator (~140KV, 200–400W) | ~$50 | Select after RPM/voltage analysis |
| Schottky diode (10A, 60V) | ~$5 | Battery bus protection |
| Belt/coupling hardware for generator shaft | ~$20 | TBD by mechanical design |
| **Phase 2 total (est., excl. labor)** | **~$1,240** | Excludes VTOL battery (Phase 1 item) |

---

## 6. Phase 3: BVLOS (Stubs)

Target: November 2026–January 2027.

### 6.1 Long-Range Telemetry (> 200 km)

E-REQ-07 requires real-time command and telemetry beyond 200 km. This rules out standard RF at any practical power level (100 km practical ceiling for RFD900x with directional antennas).

Candidate systems:

| System | Range | Latency | Cost | Weight | Notes |
|---|---|---|---|---|---|
| RFD900x (1W, 900 MHz) + high-gain tracker antenna | ~80–120 km | < 100 ms | ~$300 + tracker | Moderate | Requires tracking ground station; range still short of 200 km |
| Silvus SC4200 MIMO | ~150 km | < 50 ms | ~$2,000+ | ~300g | High cost |
| Iridium 9603 SBD | Global | 2–20 s | ~$400 + service | ~11g | Low-rate only; not suitable for real-time control |
| Iridium Edge Solar | Global | 300 ms–2 s | ~$600 + service | ~90g | Better latency; still marginal for control |
| RockBLOCK 9603 + RC link hybrid | Global telemetry + local RC backup | Mixed | ~$450 | ~20g | Good for BVLOS monitoring with local landing control |
| Starlink Mini | Global | < 50 ms | $599/mo + $599 hw | ~1 kg | Too heavy for 25 kg platform |

Decision deferred to Phase 3 design. Likely architecture: satellite link for BVLOS telemetry and command, secondary 900 MHz RF for takeoff/landing proximity ops. Confirm Turkish regulatory requirements for BVLOS satellite operation.

### 6.2 Obstacle Avoidance

E-REQ-06 requires a front-facing sensor. ArduPilot supports BendyRuler avoidance with PROXIMITY class sensors.

Candidates: RPLidar S2L (360° 2D, 32 m, USB, ~200g), NanoRadar MR82 (80 GHz, 30 m, CAN, ~65g), Intel RealSense D435 (stereo vision, USB, ~72g). Details TBD in Phase 3 design.

---

## 7. Phase 4: Payload Bay (Stubs)

Target: February–March 2027.

### 7.1 Payload Bay Power

| Parameter | Requirement |
|---|---|
| Voltage | 12V regulated (E-REQ-09) |
| Current | TBD based on payload types (surveillance cameras, cargo sensors) |
| Protection | Per-payload fuse on 12V rail |
| Source | Dedicated 12V DC-DC converter from HV bus (separate from avionics BECs) |
| Connector | Quick-release (see §7.3) |

### 7.2 Payload Data Interface

| Interface | Use |
|---|---|
| CAN bus spur | Low-rate sensor data, DroneCAN-compatible sensors |
| Ethernet | High-bandwidth: video, LiDAR point cloud. Via companion computer or Pixhawk TELEM2 → Ethernet bridge. |
| Serial (UART) | Fallback for simple payload sensors |

### 7.3 Quick-Release Connector

Arrow payload standard connector (per project requirements §9.7). Candidate: Molex Mini-Fit Jr multi-pin or custom interface combining 12V power and CAN/Ethernet in a single latching connector. Final spec in Phase 4 design, aligned with Arrow platform standard.

---

## 8. Work Packages

### WP-E01 — Battery and BEC Sizing and Procurement

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | MAD Motor spec sheet (motor efficiency vs. throttle curve), MTOW 25 kg, hover time target, CG constraints from Alperen |
| Outputs | Battery specification (capacity, C-rating, weight, dimensions, brand), BEC selections, recommended battery placement |
| Dependencies | MAD Motor spec sheet receipt (in transit); Alperen's fuselage CG envelope |
| Deliverables | Written battery spec with sizing rationale, BEC datasheets, procurement order |

### WP-E02 — HV Wiring Harness Design

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | ESC and motor physical locations (from Alperen's boom and fuselage design), battery placement, peak current from WP-E01, connector selections |
| Outputs | Harness schematic, routing diagram with cable lengths, AWG selection table, connector BOM, assembly instructions |
| Dependencies | WP-E01 complete; Alperen's boom and fuselage geometry finalized |
| Deliverables | Harness diagram (KiCad schematic or detailed hand drawing), wire and connector BOM, assembly guide |

### WP-E03 — Tail PCB Design and Fabrication

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | Here4 second cable pinout (confirmed by WP-E04), servo voltage requirement (7.8V), tail bay dimensions from Alperen, servo connector type from WP-E06 |
| Outputs | KiCad schematic and layout, Gerbers, assembled and tested PCB |
| Dependencies | **Blocked on WP-E04** (Here4 PWM validation must pass before finalizing PCB pinout); WP-E06 for servo connector type; Alperen's tail bay dimensions |
| Deliverables | KiCad project files in repo, Gerber files, 5 assembled boards, test report confirming PWM output timing and voltage |
| Note | If WP-E04 fails, pivot to Matek CAN2PWM and redesign PCB as a voltage reg and power breakout only |

### WP-E04 — Here4 CAN-to-PWM Validation

| Field | Value |
|---|---|
| Owner | Thomas (executor), Zeynep (procedure author) |
| Phase | 1 |
| Inputs | Here4 with both cables connected, breakout board from Thomas's kit, ArduPilot-configured FC (existing Quiver hardware acceptable), oscilloscope or servo response check |
| Outputs | Pass/fail: confirmed PWM output via Here4 breakout board is functional, at correct timing, and at 5V logic level |
| Dependencies | Thomas's Ankara visit (~May 12, 2026); Zeynep's test procedure document |
| Deliverables | Short test report with scope screenshots, posted to Spearhead forum |
| Note | If Thomas's kit lacks a breakout board: Alperen purchases locally in Turkey (~$600) or orders from EU supplier. Escalate to Erick. |

### WP-E05 — Avionics and Sensor Integration

| Field | Value |
|---|---|
| Owner | Erick (integration), Zeynep (ArduPilot configuration) |
| Phase | 1 |
| Inputs | Pixhawk 6C pinout, Here4 CAN pinout, altimeter selection (Ainstein vs Benewake), pitot tube model from storage, CAN bus topology from §4.9 |
| Outputs | Full avionics wiring diagram, ArduPilot parameter file (.param) covering all sensors, validated bench-tested avionics stack |
| Dependencies | WP-E04 (Here4 CAN confirmed); altimeter selection is part of this WP |
| Deliverables | Avionics wiring diagram, parameter file in repo, bench test report with GPS fix confirmation, altimeter readout, and airspeed sensor calibration |

### WP-E06 — Servo Selection and Control Surface Wiring

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | Zeynep's control surface sizing study (final torque requirements), hinge geometry from Alperen, servo bay and linkage dimensions from Alperen, wing panel connector standard |
| Outputs | Servo model selections for V-tail (×2) and flaperon (×2) with torque analysis, wing panel connector spec, servo wiring diagram |
| Dependencies | Zeynep's control surface sizing final output; Alperen's tail and wing servo bay dimensions |
| Deliverables | Servo selection with torque safety factor analysis, wing panel connector spec (connector model, pinout, drawing), servo wiring diagram |

### WP-E07 — Physical Electrical Layout and Avionics Bay Design

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | Component dimensions (Pixhawk 6C, BECs, power module, RC receiver, telemetry radio), nose bay dimensions from Alperen, battery bay dimensions, CG target from Alperen |
| Outputs | Physical layout drawing (top view and side view), cable routing plan, component mass and CG contribution |
| Dependencies | WP-E01 (battery dimensions locked); Alperen's fuselage truss design finalized |
| Deliverables | Dimensioned layout drawing (Inventor, Fusion, or dimensioned sketch), component placement table with CG contribution |

### WP-E08 — RC and Telemetry Link Selection and Integration

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | Phase 1 range requirement (VLOS, ≤ 5 km), team's transmitter hardware, GCS preference (Mission Planner or QGroundControl), ArduPilot SERIAL configuration |
| Outputs | RC system and telemetry radio selections, FC UART configuration, antenna placement plan, range test results from bench |
| Dependencies | Team transmitter and GCS confirmation from Zeynep/Alperen |
| Deliverables | Hardware selection with rationale, parameter changes, bench radio-link verification |

### WP-E09 — Heading LED System

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | Wing tip geometry from Alperen, ArduPilot NTF documentation, available GPIO pins on Pixhawk 6C |
| Outputs | LED strip type and length, wiring plan, ArduPilot NTF parameter configuration |
| Dependencies | Wing tip dimensions confirmed by Alperen |
| Deliverables | Wiring diagram, NTF parameter configuration, installed and tested (confirmed LEDs respond to arm/disarm and mode changes) |

### WP-E10 — IC Engine Electrical Integration (Phase 2)

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 2 |
| Inputs | Engine selection (DLE vs Stinger, 35cc vs 55cc), ignition module spec, throttle servo travel, starter spec if procured |
| Outputs | Phase 2 additions to harness diagram, ignition supply regulator design, throttle servo wiring, vibration isolation plan, notch filter tuning procedure |
| Dependencies | Engine selection by Alperen; Phase 1 electrical installation complete |
| Deliverables | Updated harness diagram, ignition supply schematic, Phase 2 parameter additions (.param delta from Phase 1 file) |

### WP-E11 — Generator Architecture Research and Selection (Phase 2)

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 2 |
| Inputs | Cruise power budget (100–200 W target), weight budget (~1 kg), engine shaft output and RPM at cruise, starter motor specs if Pilot RC unit is procured |
| Outputs | Generator architecture decision with supporting analysis, preliminary schematic for chosen option |
| Dependencies | Engine selection; Phase 1 flight data for actual cruise power consumption |
| Deliverables | Architecture decision document with power budget, efficiency estimate, and weight breakdown; preliminary schematic |

### WP-E12 — Long-Range Telemetry Research and Selection (Phase 3)

| Field | Value |
|---|---|
| Owner | Erick and Zeynep |
| Phase | 3 |
| Inputs | 200 km range requirement (E-REQ-07), latency tolerance for BVLOS command loop, Turkish and international BVLOS regulatory requirements, Phase 2 flight data |
| Outputs | Telemetry system selection, link budget analysis, integration and antenna plan |
| Dependencies | Phase 2 flight testing validating platform readiness for BVLOS |
| Deliverables | Link budget spreadsheet, hardware selection with rationale, wiring and parameter changes, regulatory compliance note |

---

## 9. Open Questions

| ID | Question | Owner | Priority | Gates |
|---|---|---|---|---|
| OQ-01 | Battery: exact capacity, brand, and CG position locked | Erick | High | All Phase 1 work |
| OQ-02 | Here4 CAN-to-PWM: does the PWM breakout work reliably? | Thomas / Zeynep | High | WP-E03 |
| OQ-03 | Altimeter: Ainstein US-D1 vs Benewake TF03-180? | Erick | Medium | WP-E05 |
| OQ-04 | RC system: ELRS vs Crossfire? (depends on team transmitter) | Alperen / Erick | Medium | WP-E08 |
| OQ-05 | HV kill switch: contactor + key-switch or Anderson PP75? | Erick | Medium | WP-E02 |
| OQ-06 | Pitot tube (from storage): model and interface confirmed? | Alperen | Medium | WP-E05 |
| OQ-07 | Shared or separate battery for VTOL vs IC ignition in Phase 2? | Erick | Low | Phase 2 design |
| OQ-08 | Electric starter generator mode: what is the mechanical limitation? | Erick | Low | WP-E11 |
| OQ-09 | Wing panel connector: Molex Mini-Fit Jr confirmed, or alternative? | Erick | Medium | WP-E06 |
| OQ-10 | Payload quick-release connector: align with Arrow platform standard? | Erick / Alperen | Low | Phase 4 design |
| OQ-11 | Here4 breakout board: does Thomas's kit include one? | Thomas | High | WP-E04 escalation |

---

## 10. Revision History

| Rev | Date | Author | Changes |
|---|---|---|---|
| 0.10 | May 2026 | errrks | Added STATUS callout at the top of §5 (Phase 2): the entire section is preliminary, not a design baseline. Listed the specific items requiring further research before commitment (engine selection, starter, VESC dual-role viability, generator topology, ignition battery, CHT/EGT, vibration spectrum, throttle servo). Phase 1 (§4) remains the working baseline. |
| 0.9 | May 2026 | errrks | Corrected §3 requirements traceability source: replaced "Charter §X.Y" citations with "AIP-006 §X.Y" (DAO forum proposal, March 2026). The formal project charter referenced in proposal §5 has not yet been published, so AIP-006 is the authoritative requirements source. Added clarifying note above the table. |
| 0.1 | May 2026 | errrks | Initial draft — Phase 1 specified, Phases 2–4 stubbed |
| 0.2 | May 2026 | errrks | Updated §4.1 and §4.5 with MAD V8013 PRO datasheet data: revised hover current (70A → 98–124A), battery floor (16Ah → 22Ah), added full performance table and operating zone limits |
| 0.3 | May 2026 | errrks | Updated §4.5 ESC with full AMPX 80A V2 datasheet: BEC 5V/200mA (avionics-inadequate), IP67, thermal shutdown thresholds (125°C/140°C), 2s CAN watchdog, wire specs (12AWG/800mm power, 14AWG/150mm phase). Corrected §4.6 harness AWG. Updated §4.13 with ESC thermal/CAN watchdog notes. Corrected BOM ESC unit price ($95 → $60). |
| 0.8 | May 2026 | errrks | Confirmed smart battery: removed balance tap, PM07, and current sensor from design. Updated §4.1 (smart battery note), §4.2 (distribution), §4.3 (kill switch scope), §4.6 (harness), §4.7 (FC power), §4.13 (battery monitoring via DroneCAN BatteryInfo). BOM updated — nose GPS and RC/telemetry TBD pending team decisions. |
| 0.7 | May 2026 | errrks | Corrected §4.1 hover analysis: 12S voltage range closely matches datasheet 48V conditions; motor operates in continuous zone for most of flight discharge curve, short-term zone only near low-battery threshold. Corrected capacity sizing to average draw. Separated §4.13 ESC thermal (junction temp, 125/140°C) from motor thermal (winding current zones, 24A/78A, motor datasheet source). Added motor temp monitoring proxy approach. |
| 0.6 | May 2026 | errrks | Removed motor tilt (decision reversed). Renamed V-tail surfaces to ruddervators throughout. Restructured servo power: 7.8V regulation split — flaperon BEC in CG bay, ruddervator regulation on tail PCB from HV spur (no LV servo run through fuselage). Updated tail PCB spec to HV input with on-board switching regulators. Added smart battery decision point in §4.1. |
| 0.5 | May 2026 | errrks | Updated §4.5 propeller with full FLUXER PRO 26×7.8 MATT datasheet: 69g/blade, Ø10mm bore, Ø20mm M3×4 hub, 1900–4000 RPM optimum, 18 kgf thrust limit, prop test data cross-validated against motor datasheet. Corrected BOM prop price ($90 → $193.90/pair). Phase 1 grand total updated ($4,668 → $5,081). |
| 0.4 | May 2026 | errrks | Expanded §5 Phase 2 from stubs to detailed architecture: fuel mix and storage (§5.1), throttle and ArduPilot ICE parameters (§5.2), choke actuation and cold-start sequence (§5.3), ignition and kill wiring (§5.4), electric starter circuit (§5.5), RPM sensor options and ArduPilot notch filter setup (§5.6), engine health monitoring CHT/EGT (§5.7), VESC dual-role starter-generator architecture with battery top-up design (§5.8), vibration isolation procedure (§5.9), Phase 2 BOM (§5.10). |
