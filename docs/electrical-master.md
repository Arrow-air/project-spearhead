# Spearhead Electrical Master Document

**Doc ID:** SPH-E-001
**Revision:** 0.16
**Author:** errrks
**Date:** June 2026
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
   - 4.6 Wiring Harness and Connection Plan
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
                              │  XT150 (pack socket, ground safing disconnect)
                  [HV Kill (PT1 implementation TBD, OQ-05)]
                              │  8 AWG
                        ┌─────┴──────┐
                        │ FCHUB-12S  │  44.4V nominal
                        │    PDB     │  battery pads = FC power tap point
                        └┬───┬───┬──┬┘
                        ESC1 ESC2 ESC3 ESC4  (all DroneCAN)
                          │   │   │   │
                          M1  M2  M3  M4   (FR / RL / FL / RR)

   FCHUB battery pads ─ 16 AWG XT60 pigtail ─ [PM02D HV]
                                                5.2V/3A + I2C → 6-pin → Pixhawk POWER1

                  FCHUB-12S PDB rails
                        │
          ┌─────────────┼──────────────────┐
          │             │                  │
   [5V/5A BEC]  [HV spur → Tail PM12S-3]  [12V/4A BEC]
  peripheral 5V   (Vx → ruddervators,     payload (Ph 4)
  (CAN/RC/telem)   5V → adapter+Here4)    + contactor coil
                                           if D2 picks one)

  Flaperon servo power: forward PM12S-3 Vx rail (HV input
  at the FCHUB pads), not shown above. See §4.4.
  Tail PM12S-3 HV input: 18 AWG branch off the ESC-lead extension
  joint at the boom root (not a FCHUB rail), behind the HV kill. See §4.10.
          │
    [LV Kill Switch]
   (FC + periph 5V)
          │
     ┌────┴──┬────────┐
  Pixhawk   GPS   Sensors/RC/
    6X    (shelf) Telemetry
 (POWER1 = power module, POWER2 unused, single FC supply)
```

### 2.2 Communications Topology

```
Pixhawk 6X
│
├── CAN1 ─── [120Ω] ── CAN splitter/hub
│        ├── ESC1..ESC4 (nodes 1–4, short boom stubs)
│        ├── GPS #1 (fwd shelf, node 10, primary, model TBD)
│        └── Here4 (tail, node 11) ── CAN-PWM adapter (node 12) ── [120Ω]
│
├── CAN2 ─── [reserved / Phase 3 expansion]
│
├── UART1 (TELEM1) ─── Ground telemetry radio (MAVLink 2.0)
├── UART2 (TELEM2) ─── [reserved / Phase 3 long-range link]
├── UART3 (GPS1)   ─── [unused, GPS on CAN]
├── UART4          ─── Radar/LiDAR altimeter
├── UART5          ─── RC receiver (CRSF or SBUS)
│
├── I2C   ─── 
├── ADC   ─── Pitot static pressure transducer
│
├── PWM OUT (IO, outputs 1–8):
│     Output 1: Aileron, right flaperon (physical PWM from the 6X IO)
│     Output 2: Aileron, left flaperon  (physical PWM from the 6X IO)
│     Output 3: RuddervatorLeft  (SERVO_FUNCTION 79, via DroneCAN to tail adapter, node 12)
│     Output 4: RuddervatorRight (SERVO_FUNCTION 80, via DroneCAN to tail adapter, node 12)
│     Outputs 5–8: VTOL motor commands via DroneCAN ESC nodes 1–4 (no physical PWM)
│
└── USB-C ─── Ground configuration / parameter upload
```

VTOL motor commands go over DroneCAN via ESC_BM bitmask. Physical PWM outputs 5–8 are unused in Phase 1 unless DroneCAN falls back to PWM.

CAN bus length nose to tail: approximately 2.5 m. No signal integrity issue at 1 Mbps for this length.

### 2.3 Physical Location Overview

| Zone | Contents |
|---|---|
| Nose (slimmed, v5) | No avionics requirement. Pitot mast and external antenna placement only if needed |
| Forward battery bay | ProFuse 12S 16Ah 60C pack (3,822g, 186×76×160mm), XT150 pack connector (ground safing disconnect), HV kill (implementation TBD, OQ-05). Moved forward for CG per the v5 layout |
| Avionics shelf (above battery, top access) | Pixhawk 6X, PM02D HV power module (FC POWER1, single supply), primary GPS on the shelf cover (sky view), RC receiver, telemetry radio, LV kill switch |
| CG bay / boom roots | Matek FCHUB-12S PDB (HV distribution + 5V/12V BECs, battery pads carry the FC power tap), CAN splitter/hub, forward PM12S-3 (flaperon Vx rail) |
| Boom roots / tips (×4) | ESC externally mounted on boom near motor (within 150mm of motor centerline) |
| Motors (×4) | MAD V8013 PRO IPE 150KV, 3–5° forward tilt relative to aircraft body (motors vertical relative to gravity when aircraft is in cruise attitude) |
| Tail group (boom-mounted) | Here4 (secondary GPS), Matek CAN-L4-PWM CAN-to-PWM adapter, PM12S-3 (Vx servo rail + fixed 5V for adapter and Here4), ruddervator servo connectors |
| Wing panels (×2) | Flaperon servos (5V for PT1), wing panel connector (signal + 5V + GND, 5V from the peripheral rail) |
| Airframe exterior | WS2812 heading LED strips (wing tips, nose, tail) |

Zones follow Alperen's fuselage conceptual v5 (June 11): longer thinner nose, forward battery for CG, avionics shelf above the battery with top access, pusher between the fuselage and the boom-carried tail. The tail harness (CAN trunk + HV spur) therefore routes along a tail boom, not through a fuselage run. Shelf and bay dimensions are open, see `docs/open-questions.md`.

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

**Battery capacity sizing:** the QuadPlane hovers only during VTOL takeoff, landing, and transition and flies wing-borne the rest of the mission, so sustained hover is not the energy driver. At ~85A average hover draw a 16 Ah pack gives roughly 11 minutes of hover before reserve, well beyond the per-flight hover budget for Phase 1 testing and transitions. The bay is sized to **16,000 mAh** (12S) per Alperen's mission-energy rerun (May 2026), trading the earlier conservative 22 Ah sustained-hover figure for about 1.7 kg of weight. Confirm against Alperen's published mission-energy note.

**Peak current** (safety-margin thrust at 75% throttle, 36.6A/motor × 4):
- 4 × 36.6A = **146A** (short-term zone, < 30s)
- At 16,000 mAh: 146 / 16 = **9.1C** (within high-discharge LiPo ratings)

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
| Capacity | **16,000 mAh** (12S, per Alperen's mission-energy rerun) |
| C-rating | 60C (ProFuse spec), far above the ≥15C floor. The 320A peak draw is trivial for the pack |
| Connector | XT150 socket on the pack (ProFuse-fitted). Harness side: XT150 male, 8 AWG to the FCHUB battery pads (4 AWG does not fit XT150 solder cups, see §4.6). Plain XT150 is not anti-spark, see the kill interaction note in §4.3. No Y-splitter: the FC power module taps the FCHUB pads |
| Form factor | **Selected June 11: ProFuse Super Nano 12S 16,000 mAh 60C**, single pack, 186×76×160mm, robotsepeti (TR, in stock, 41,322 TL incl. VAT). Forward battery bay per the v5 fuselage, the 160mm height drives the bay design. Lock position with Alperen before the fuselage is built |
| Weight | 3,822g (ProFuse datasheet), largest single CG impact item |

**Dumb LiPo pack (no BMS or CAN).** A standard high-discharge 12S LiPo is used, not a smart battery. Alperen found the smart-battery option added roughly 1 kg, not worth it for the end product. Implications:
- No CAN battery telemetry and no per-cell data. Vehicle current comes from the FCHUB-12S 440A sensor and pack voltage from the power module or a FC voltage sense (see §4.2 and §4.13)
- Balance leads are used for charging only, not wired into the aircraft. In flight the pack connects to the HV path by main leads only
- Set conservative pack-level low-voltage thresholds (no per-cell cutoff). Carry a plug-in cell checker for ground checks
- Pack selected June 11 (WP-E01 selection closed): ProFuse Super Nano 12S 16Ah 60C, XT150 socket

---

### 4.2 Power Distribution

| Parameter | Value |
|---|---|
| Board | Matek FCHUB-12S **V2**, off-the-shelf PDB (no custom board for V1) |
| Architecture | Star: battery to the FCHUB-12S, individual leads to each ESC. The battery input pads also carry the FC power module tap (16 AWG XT60 pigtail), replacing the earlier battery Y-splitter |
| Input | 8–60V (3–12S), battery main leads only (no balance harness) |
| Per-channel rating | 70A continuous, 110A burst per ESC pad. Hover is ~24–31A per motor, peak ~80A per motor (capped by the 80A ESC), so the margin is ample |
| ESC connection | 4× ESC power pads. Only the power pads are used. The signal and telemetry pads stay unused because the ESCs run DroneCAN (see §4.9) |
| On-board BECs | 5V/5A (max 6A), 12V/4A (max 5A), 3.3V/0.5A. Rail assignment in §4.4 |
| Current sensor | On-board 440A analog sensor used for vehicle current (3.3V ADC, scale 75, to a FC ADC pin). Pack voltage from the V2 built-in 1K:20K divider (V pad set to VBat). See §4.13 |
| V2 breakout | JST-SH 8-pin: V (selectable 12V or VBat via the 1K:20K divider), G, Curr, TLM (BLHeli32, unused), S1–S4 (ESC signal, unused with DroneCAN). Ships with the SH1.0 8-pin cable |
| Switchable 12V | 12VSW pad, on/off via FC PINIO (GPIO level), 2A constant. Candidate supply for the EV200 contactor coil (§4.3) or a switched LED feed |
| Per-branch fusing | 100A fast-blow per ESC branch, optional for Prototype 1 (assess in WP-E02) |
| Physical | 55 × 50 × 6 mm, 30.5 × 30.5 mm Φ3 mounting, 21g. Ships with a low-ESR bulk capacitor for ESC switching transients |

The FCHUB-12S replaces the earlier custom PCB and bus bar option, in line with the V1 goal of minimizing custom boards. Holybro PM07 is not needed since the FCHUB current sensor handles vehicle-current monitoring. The 70A per-channel continuous rating comfortably exceeds the AMPX 80A ESC ceiling, and the 4×110A burst rating covers the 320A peak.

---

### 4.3 Kill Switches

| Switch | Covers | Does NOT cover | Placement | Current rating |
|---|---|---|---|---|
| HV Kill | 12S battery → FCHUB PDB (including the FC power module pad tap), tail BEC HV spur | — | Exterior, reachable within 3 s, plus RC-switchable per the June 11 criterion | ≥ 300A peak capability |
| LV Kill | FC power (power module 5.2V to POWER1), peripheral 5V rail (GPS, RC, telemetry, flaperon servos) | Tail PM12S-3: ruddervator servos + adapter + Here4 (on HV spur) | Exterior, separate from HV | ≥ 10A |

The tail HV→5V BEC powers the ruddervator servos, the CAN-PWM adapter, and the Here4, all from the HV bus spur. The whole tail control path therefore has no separate LV kill and de-energizes only when the HV kill is thrown, preserving ruddervator authority whenever the aircraft is armed. This is intentional, and it resolves the earlier concern about the adapter losing power on an LV kill, since it no longer draws from the CAN 5V rail.

**PT1 implementation (open, criterion added June 11):** AIP-006 §9.3 SHALL requires independent HV and LV kill switches, and those drive the design above. The June 1 review noted the contactor approach came from Quiver and final-product assumptions and may be heavier than PT1 needs. The June 11 review added a criterion: the HV kill must be switchable from the RC link, not only from a physical disconnect at the airframe. That rules out a pure manual-disconnect deviation and narrows PT1 to three candidates (decide in OQ-05 before the PT1 order ships):

| Candidate | How it switches | Mass | Notes |
|---|---|---|---|
| TE Kilovac EV200AAANA contactor | 12V economizer coil (1.7W hold ≈ 0.14A) from the FCHUB V2 12VSW pad, switched directly by FC PINIO. No driver hardware. RC kill via RC_OPTION relay mapping, GCS kill via MAVLink relay | ~450g | ~€149. Proven interrupt rating (2000A at 320VDC). Closing onto the discharged ESC input caps at 50V is within rating, no precharge stage needed |
| MOSFET anti-spark switch (e-skate class, 200A cont / 400A peak, up to 20S) | Enable line interfaced to an RC switch or FC GPIO | ~100–150g | ~$80. Fails short (a failed kill no longer kills), no flight heritage, enable interface is an adaptation. Bench evaluation required before trusting it in the HV path |
| Software E-stop + ground safing only | ArduPilot motor emergency stop on an RC switch cuts CAN motor commands. AMPX throttle-loss protection backstops at 2 s. XT150 unplug for ground safing | 0g | Not a hardware kill. A recorded §9.3 deviation for PT1 only, full hardware kills return for the final product |

The software E-stop layer costs nothing and gets configured regardless of which hardware path wins. The earlier WP-E02 options table (contactor + key switch, Anderson PP75 pull disconnect, E-stop pushbutton) is superseded: none of those are RC switchable, and the PP75 was undersized anyway at 75A-class continuous against a 98–124A hover draw.

**Connector interaction (June 11):** the selected ProFuse pack carries a plain XT150 socket, not an anti-spark connector. With an inline hardware kill (contactor or MOSFET switch) the battery is connected with the kill open, so there is no mating spark and the kill element absorbs the ESC cap inrush. That is one more argument for a hardware kill. If the software-only path wins, add an AS150U anti-spark adapter pigtail at the pack or accept the mating spark at 12S.

LV kill switch: illuminated rocker or key switch, 10A @ 12V. Interrupts the FC power (power module output) and the peripheral 5V rail (GPS, RC, telemetry, flaperon servos). It does not cut the tail PM12S-3, so the tail control path stays powered on the HV spur until the HV kill is thrown.

---

### 4.4 Voltage Regulation

| Rail | Voltage | Continuous Current | Loads | BEC type |
|---|---|---|---|---|
| FC power (POWER1) | 5.2V | 3A max | Pixhawk 6X (single supply, POWER2 unused) | PM02D HV (I2C, in the 6X set), fed from a 16 AWG XT60 pigtail at the FCHUB battery pads |
| Peripheral 5V | 5.0–5.2V | 4–6A | Nose/shelf GPS (CAN bus), RC receiver, telemetry radio, sensors | FCHUB 5V, separate from the FC inputs. Tail devices (adapter, Here4) run on the tail PM12S-3 instead |
| LED 5V | 5V | up to ~3A | WS2812B strips (wing tips, nose, tail) | Separate feed, high transient draw, kept off the FC rail |
| Flaperon servo | 5.25V (PT1), 8V (final) | 2–4A PT1, 15A available | 2× flaperon servos on wing panels | Forward PM12S-3 Vx rail via the wing connectors. PT1 test servos, see §4.11. |
| Ruddervator servo | 5.25V (PT1), 8V (final) | 2–4A PT1, 15A available | 2× ruddervator servos in tail bay | Tail PM12S-3 Vx rail, fed from the HV bus spur. Its fixed 5V powers the adapter and Here4. No long LV run through the fuselage. |
| Payload | 12V | TBD | Phase 4 only — not needed for Phase 1 | — |

FC power: the Pixhawk runs from a single source, the power module on POWER1. POWER2 is unused for V1, so there is no hot-standby. This is a single point of failure on FC power, accepted to keep V1 simple. Add a POWER2 backup BEC later if FC power redundancy is wanted before extended flights.

**FC supply (Pixhawk 6X selected, June 11):** the FC runs from the PM02D HV (I2C, 2–12S, included in the 6X standard set) on POWER1, fed by a 16 AWG XT60 pigtail soldered at the FCHUB-12S battery input pads. The pads are the same electrical node a battery Y-splitter would create, with one fewer connector pair in the 320A path, so the earlier AS150U/XT150 Y-splitter is dropped (it remains the fallback if the pad tap proves mechanically awkward). For the record: the PM02D is supported on the 5X/6X series only, the abandoned 6C path would have needed the analog PM02 V3 (verified against Holybro, June 2026). The module carries only the FC draw.

**FCHUB-12S on-board BECs:** the FCHUB provides 5V/5A (max 6A), 12V/4A, and 3.3V/0.5A (verified against the Matek page, June 2026). With the FC on its own power module, use the FCHUB 5V for the peripheral 5V rail (CAN GPS, adapter, RC, telemetry), and the FCHUB 12V for the Phase 4 payload and the HV kill contactor coil if an EV200 or GX11 is chosen. Keep the WS2812 LED strips on a separate 5V feed so their transients stay off the peripheral rail.

BEC selection (updated late June 11): **2× Matek PM12S-3** (Aykut Havacılık TR, 6,900 TL each, specs verified June 2026). Per unit: 9–55V input with TVS protection, Vx rail selectable 5.25V / 6V / 8V at 15A continuous (25A peak), fixed 5.2V/4A, fixed 12V/4A, built-in 1K:20K divider, 53g. Roles: the **tail unit** Vx feeds the ruddervator servo rail (5.25V for PT1, 8V for the final HV servos) and its fixed 5V feeds the adapter logic and the Here4. The **forward unit** Vx feeds the flaperon servo rail through the wing connectors, and its 12V is a second candidate supply for the contactor coil. The 15A Vx rating removes servo current headroom questions, and the later 8V change is a pad setting, not new hardware. Supersedes the BEC12S-PRO plan and the earlier candidates (Matek MBEC6S, Castle BEC Pro, Holybro PM07).

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
CAN_D1_UC_SRV_BM = 0x000C   ; servo outputs 3–4 (ruddervators) to the tail adapter, node 12
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

### 4.6 Wiring Harness and Connection Plan

#### Wire Sizing

Current ratings below assume silicone-insulated stranded copper in open-air UAV duty (not bundled in conduit).

| Segment | AWG | Est. Length | Est. Current | Connector |
|---|---|---|---|---|
| Battery main leads | **8 AWG** | 0.3–0.6 m (forward bay to FCHUB) | ~98–124A hover, 320A peak | XT150 male at the pack socket (cups take 8 AWG, not 4 AWG), leads land on the FCHUB battery pads. FC power module taps the same pads via a 16 AWG XT60 pigtail |
| HV bus to each ESC (×4) | **12 AWG** | ~1.4 m run (rear booms 20–30 cm longer) | ~25–31A hover, 80A rated peak | XT90S or solder direct to ESC power leads |
| ESC power leads (factory, ESC to FCHUB) | **12 AWG** | **800mm** (factory) | 80A rated | Bare tinned. Add XT90S or solder to FCHUB pads |
| ESC phase wires to motor (factory) | **14 AWG** | **150mm** (factory) | 80A rated (short duty only) | Bare tinned — 5mm bullet at motor side |
| BEC input from HV bus | 16 AWG | 0.2–0.3 m | 5A max | XT30 |

**Wire notes:**

The AMPX 80A ESC ships with 12AWG/800mm power leads and 14AWG/150mm phase leads. These are MAD's designed wire lengths for their thermal model at rated current.

- **Phase leads (14AWG/150mm): do not extend** without upsizing to 12AWG or heavier. The 150mm is short by design — heat dissipates through the motor housing. A longer 14AWG run at 80A will overheat.
- **Power leads (12AWG/800mm):** adequate for Spearhead's hover current (~25–31A/motor). At rated 80A the wires run warm but duty cycle is short (10–30s per motor zone limits). The boom run exceeds 800mm (confirmed June 11), so each ESC power lead is extended at the boom-fuselage junction with an XT90S splice. Extend with **12AWG throughout** to match the factory lead. MAD rates that same 12AWG factory lead at the full 80A, and the extension carries the same current in series while Spearhead never sustains more than ~31A/motor in hover (the motor winding limit caps continuous draw at 24A, §4.5), so no upsize is warranted. Matching 12-to-12 also keeps the splice a clean same-gauge joint instead of a 12-to-10 step.
- **Main battery leads:** 8AWG silicone with XT150. The earlier 4AWG spec does not fit XT150 solder cups, and on this 0.3–0.6 m run the connector (~130–150A continuous class), not the wire (~150A for 8AWG silicone in open air), is the thermal limit. Hover at 98–124A runs the lead warm, the 320A peak lasts seconds. Monitor connector and lead temperature during hover tests. The FC power module leg is small-gauge (FC supply only), an XT60 pigtail soldered at the FCHUB battery pads.

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
| HV spur to tail | 18 AWG HV spur broken out at the tail-boom ESC-lead extension joint (boom root) to the tail BEC (5V is regulated at the tail, not run from the nose) |
| Signal logic level | 5V (from the CAN-PWM adapter) |
| Power logic | 5V (local tail BEC, shared ground) |

#### Wing Panel Connector

The detachable wing panels each carry one flaperon servo. The connector at the wing root must break out signal, servo power, and ground.

| Signal | Notes |
|---|---|
| Servo signal (×1 per wing) | 5V logic PWM from nose FC |
| 5V servo power | From the peripheral 5V rail via fuselage harness |
| Ground | Shared with avionics ground |

Connector candidate: Molex Mini-Fit Jr 6-pin (known from Quiver, rated 9A per pin, latching). TBD in WP-E06.

#### Routing and EMI

Servo PWM wires run near ESC DC power in places. Alperen flagged a ~1 m parallel run between a main-wing servo PWM line and ESC power. Guidance for the discrete harness:
- Separate signal and power runs where practical. A few centimeters of gap over a 1 m parallel run is usually enough at these currents.
- Where a signal line must cross power, cross at right angles rather than running parallel.
- Use twisted pair for PWM (signal with its own return) to reject common-mode pickup. CAN is already a twisted pair (see the CAN wiring spec above).
- Keep PWM and CAN away from the motor phase wires, which carry the switching currents and radiate the most. Phase wires stay short and at the motor by design.
- For a long run that must sit near power, use shielded cable with the shield grounded at one end only, the FC end. The same single-end grounding applies to the Phase 2 CDI shielding (§5.4).
- Quiver avoids this class of noise mostly through multi-layer PCB routing with ground and power planes. Spearhead's discrete wiring does not have that, so the separation and twisting above carry the load.
- The tail PWM is generated locally at the adapter (§4.10), so only short tail-bay servo leads run as PWM, not a 2.5 m run through the ESC field.
- Flaperon PWM runs physically from the 6X IO outputs to the wing roots (~1–1.5 m per side). These are exactly the runs the separation, twisting, and crossing rules above exist for. A wing CAN-PWM node was considered on June 11 and dropped the same day with the 6X selection.

#### PT1 Harness Connection Table

Every PT1 cable end to end, per the June 11 request for the general layout. Lengths are estimates against the v5 fuselage and get confirmed once Alperen publishes shelf and bay dimensions. EMI rules from the previous subsection apply to H8, H10–H13.

| # | From | To | Function | Wire | Connector(s) | Est. length | Notes |
|---|---|---|---|---|---|---|---|
| H1 | Battery +/− | FCHUB-12S battery pads | 12S HV main | 8 AWG | XT150 male at the pack socket, solder at pads | 0.3–0.6 m | HV kill (OQ-05) sits in this run. XT150 is not anti-spark: connect with the kill open |
| H2 | FCHUB battery pads | PM02D HV input | FC power leg | 16 AWG | XT60 pigtail | 0.2 m | Pad tap replaces the Y-splitter |
| H3 | FCHUB ESC pads ×4 | ESC power leads | HV per motor | 12 AWG (factory 800mm + 12 AWG extension) | Solder at pads, XT90S splice at the boom-fuselage junction (boom exceeds the 800mm factory lead, confirmed June 11) | per boom geometry | Rear booms run longer, exact extension length pending Alperen's boom dimensions. On the tail-carrying boom this extension joint also breaks out the 18 AWG tail HV feed (H5) |
| H4 | ESC ×4 | Motor phase leads | 3-phase | 14 AWG (factory 150mm) | 5mm bullets | 150mm | Do not extend (wire notes above) |
| H5 | ESC-lead extension joint (boom root) | Tail PM12S-3 input | Tail HV spur | 18 AWG | Branch broken out at the ESC-lead extension joint (XT90S or soldered lap, per H3), XT30 at the PM12S-3 | ~1.5–2 m (boom-root joint to the tail group) | Breaks out at the same boom-fuselage junction where the tail-boom ESC lead is extended (H3), so no separate splice into a healthy conductor. Still behind the HV kill (downstream of the FCHUB ESC pad). Runs aft along the tail boom, path TBC with Alperen (inside-tube routing pending structural check) |
| H6 | PM02D HV 6-pin | Pixhawk POWER1 | 5.2V + I2C V/I sense | Module cable | CLIK-Mate / 6-pin | 0.3 m | |
| H7 | FCHUB 5V BEC | LV kill switch, then peripheral 5V rail | RC RX + telemetry power | 20 AWG | Soldered bus / servo leads | Shelf-local | LED strips get their own feed (H14) |
| H8 | Pixhawk 6X IO PWM 1–2 | Wing root connectors | Flaperon PWM | 26 AWG twisted with its GND | Molex Mini-Fit Jr 6-pin | ~1–1.5 m per side | EMI rules apply: separation from ESC power, right-angle crossings, twisted with its own return |
| H9 | Forward PM12S-3 Vx out | Wing root connectors | Flaperon servo power | 20 AWG | Same Mini-Fit Jr | ~1–1.5 m per side | 5.25V PT1, 8V final. The forward PM12S-3 HV input taps the FCHUB pads |
| H10 | Pixhawk CAN1 | CAN hub (CG bay) | DroneCAN trunk | 24 AWG twisted pair | JST-GH 4-pin | 0.3–0.5 m | 120Ω at the FC end |
| H11 | CAN hub | ESC1–4 CAN leads | ESC stubs | Factory 1 m signal leads | JST-GH / solder | Trim toward <0.3 m where layout allows | Daisy-chain along booms instead if stubs misbehave (§4.9) |
| H12 | CAN hub | Primary GPS (shelf) | CAN branch | 24 AWG twisted pair | JST-GH | 0.3–0.5 m | Node 10 |
| H13 | CAN hub | Here4, then CAN-L4-PWM (tail) | CAN tail trunk | 24 AWG twisted pair | JST-GH daisy-chain | ~1.5–2 m | 120Ω at the last tail node. Runs alongside H5 on the boom, stays twisted |
| H14 | FCHUB 5V (separate feed) | WS2812 strips | LED power, data from FC GPIO | 22 AWG | JST | Per airframe | Transients stay off the peripheral rail |
| H15 | Tail PM12S-3 outputs | Vx → adapter servo V+ rail + ruddervators. Fixed 5V → adapter logic + Here4 | Tail servo and logic rails | 20 AWG | Servo leads / JST-GH | <0.3 m | Keep the two feeds separate, Vx goes to 8V in the final design (§4.10) |
| H16 | CAN-L4-PWM PWM 1–2 | Ruddervator servos | PWM, 5V logic | Servo leads | JR 3-pin | <0.3 m | Generated locally, no long PWM run |
| H17 | Pixhawk TELEM1 | Telemetry radio | MAVLink UART | JST-GH 6-pin | JST-GH | 0.2–0.3 m | Radio model pending OQ-04 / parts list D3 |
| H18 | Pixhawk UART5 | RC receiver | CRSF or SBUS | JST-GH | JST-GH | 0.2 m | Receiver model pending OQ-04 / parts list D3 |
| H19 | Pixhawk UART4 | US-D1 altimeter | UART | JST-GH to US-D1 lead | JST-GH | 0.3–0.5 m | US-D1 also offers CAN if the UART budget tightens |
| H20 | FCHUB V/I sense | Pixhawk 6X AD&IO port ADC pins | Vehicle current + voltage (BATT_MONITOR 4) | 26 AWG | JST | Shelf-local | V pad set to VBat (1K:20K divider, multiplier ≈ 21) |
| H21 | Pitot transducer | Pixhawk I2C or ADC | Airspeed | Per model | Per model | Shelf-local | OQ-06 |

---

### 4.7 Flight Controller

| Parameter | Value |
|---|---|
| Model | **Pixhawk 6X (Holybro Standard Set v2A, selected June 11).** Set includes the FC module, Standard Baseboard v2A, PM02D HV power module, cable set, and an M9N GPS. Selection rationale: the 6X IO PWM pins drive the wing flaperons directly, saving the second CAN-PWM node and shelf space |
| Firmware | ArduPlane 4.x with QuadPlane |
| IMU | Triple redundant, on a vibration-isolated and damped sensor board |
| Barometer | Dual internal + external baro on CAN (GPS units) |
| Magnetometer | External compass on CAN (GPS units, internal mag disabled) |
| CAN ports | FMU-CAN1 (active), FMU-CAN2 (reserved) |
| UART ports | 6× (TELEM1, TELEM2, GPS1, GPS2, SERIAL4, SERIAL5) |
| PWM outputs | 16× (8× FMU direct, 8× via IO processor) |
| Power input | POWER1 from the PM02D HV (I2C, included in the set). POWER2 unused (single FC supply for V1) |
| Vibration isolation | Internal gel-foam isolators in chassis |
| Location | Nose avionics bay, oriented with X-axis forward, mounted to isolator plate |
| USB-C | Ground access for configuration and parameter upload |

Sourcing (June 11): Aykut Havacılık lists the 6X Standard Set v2A at 38,220 TL (~$955) but is out of stock. Holybro direct is $321, roughly a third of the TR price, so order from Holybro EU unless Aykut restocks within the week. The included M9N GPS is serial (not DroneCAN) and is kept as a bench and backup unit, it is not part of the dual DroneCAN GPS architecture (§4.8). For the record: the PM02D is supported on the 5X/6X series only, and earlier PM07/PM06 references are superseded (§4.2 FCHUB-12S handles distribution and vehicle current).

---

### 4.8 GPS and Navigation Sensors

#### GPS (×2 redundant, AIP-006 §9.4)

Two independent DroneCAN GPS units. The tail unit is a CubePilot Here4 (carried over from the original tail electronics plan, now GPS only). The nose unit is the primary and its specific model is still open (see OQ-12). The M9N included in the 6X standard set is serial, not DroneCAN, and stays a bench and backup unit outside this architecture. Candidates are u-blox F9P, F9P Neo, or salvaged F9P hardware (RTK), any with internal baro and compass.

| Parameter | Primary (nose) | Secondary (tail) |
|---|---|---|
| Model | DroneCAN GPS, model TBD (OQ-12) | CubePilot Here4 |
| Interface | DroneCAN | DroneCAN |
| Internal sensors | Baro and compass (depends on module) | Barometer, IMU, magnetometer |
| RTK capable | Depends on module | Yes |
| GNSS constellations | GPS, GLONASS, Galileo, BeiDou (target) | GPS, GLONASS, Galileo, BeiDou |
| Supply voltage | 5V (CAN bus) | 5V (CAN bus) |
| CAN node ID | 10 | 11 |

ArduPilot parameters:
```
GPS_TYPE = 9    ; DroneCAN
GPS_TYPE2 = 9
GPS_GNSS_MODE = 0     ; auto
GPS_GNSS_MODE2 = 0
GPS_AUTO_SWITCH = 1   ; auto-switch to best GPS
```

**GPS yaw / heading:** moving baseline GPS yaw (dual-antenna heading) requires two matched RTK units configured as base plus rover. A mixed pair (different nose module plus tail Here4) does not support moving baseline, so Phase 1 yaw comes from the compasses on the GPS units rather than from GPS. If GPS yaw is wanted later, select two matched RTK units (e.g. two Here4, or two u-blox F9P-based units), tracked in OQ-12. The ~2 m nose-to-tail separation still gives spatial diversity for the redundant fix.

#### Airspeed Sensor

| Parameter | Value |
|---|---|
| Hardware | Pitot tube from storage (model TBD — Alperen to confirm) |
| Interface | I2C or analog ADC depending on pressure transducer |
| Location | Nose or boom, away from prop wash and wing downwash |
| ArduPilot | ARSPD_TYPE = 1 (analog) or 7 (I2C DLVR) depending on transducer |

If the stored pitot's transducer is unknown or incompatible, replacement option: Matek AP_ANALOG_AIRSPEED or Holybro airspeed sensor (I2C DLVR, ~$25).

#### Radar Altimeter — Ainstein US-D1 (selected June 11, OQ-03 closed)

| Model | Type | Range | Update rate | Interface | Supply | Weight | Cost |
|---|---|---|---|---|---|---|---|
| **Ainstein US-D1 (selected)** | Radar | 0.5–50 m | 100 Hz | UART or CAN | 5–13V (5.5V rec.) | 110g | Quote (~$280 est.) |
| Benewake TF03-180 (not selected) | LiDAR | 0.1–180 m | 100 Hz | UART | 5V | 76g | ~$250 |

Selection rationale: radar tolerates the prop wash dust kicked up during VTOL landing on dirt fields, where LiDAR returns degrade. Specs corrected against Ainstein's current page (June 2026): the US-D1 is 110g and 100 Hz, not the 56g / 12 Hz carried in earlier revisions, so the weight argument for LiDAR is gone but the robustness argument stands. Sourcing: a unit is believed to be on hand in Turkey already (June 11). Confirm its location and interface variant (UART vs CAN) before buying anything. Only if it does not turn up: no public price, quote via Ainstein or a distributor, and US sourcing plus Turkish customs would make it the longest lead item on the list (parts list line 5).

ArduPilot integration: UART on SERIAL4 with the rangefinder protocol and RNGFND1_TYPE set to the USD1 serial driver, or the DroneCAN variant if the UART budget tightens.

---

### 4.9 CAN Bus Architecture

CAN1 bus node map:

| Node ID | Device | Type | Notes |
|---|---|---|---|
| 1 | ESC1 — Motor 1 (FR) | AMPX 80A | DroneCAN ESCStatus + RawCommand |
| 2 | ESC2 — Motor 2 (RL) | AMPX 80A | |
| 3 | ESC3 — Motor 3 (FL) | AMPX 80A | |
| 4 | ESC4 — Motor 4 (RR) | AMPX 80A | |
| 10 | GPS #1 (nose) | GPS + baro + compass | Primary GPS (DroneCAN, model TBD, OQ-12) |
| 11 | Here4 (tail) | GPS + baro + compass | Secondary GPS |
| 12 | CAN-PWM adapter (tail) | Matek CAN-L4-PWM | DroneCAN→PWM for ruddervators |

**Topology:** CAN1 leaves the Pixhawk and runs to a passive CAN splitter/hub near the CG, which fans out short stubs to the four boom ESCs. The trunk continues from the splitter to the tail (Here4 and adapter) and branches to the nose GPS. This keeps the CAN fan-out off any custom board for V1. Keep each splitter stub short (ideally under 0.3 m at 1 Mbps), and confirm whether the splitter has a built-in 120Ω terminator. If the boom stubs show marginal signal integrity, daisy-chain the ESCs along the booms instead.

Bus termination: 120Ω at the Pixhawk CAN port (built-in or terminator plug) and 120Ω at the last physical node at the tail (whichever of the Here4 or CAN-PWM adapter sits at the bus end). The CAN splitter sits mid-bus and is not terminated. No termination on the other intermediate nodes.

CAN2 is reserved. Future uses: generator controller (Phase 2), companion computer CAN bridge, additional sensors.

ESC node IDs are assigned via DroneCAN UI Tool (pre-flight on bench, not in the field). Set before boom installation.

---

### 4.10 Tail Electronics

**Problem:** Running PWM servo signal cables 2.5 m from the nose avionics bay to the tail surfaces through an environment with 4× 80A ESC switching noise causes signal degradation and potential interference.

**Solution:** CAN to the tail, PWM generated locally at the tail. A standalone DroneCAN-to-PWM adapter at the tail converts the ruddervator servo commands (broadcast over CAN) into local PWM. The tail Here4 is a GPS node only. A local Matek PM12S-3 powers everything: its Vx rail feeds the ruddervator servos and its fixed 5V feeds the adapter and the Here4. No custom PCB, and no long LV run through the fuselage.

> **Supersedes:** earlier revisions used the tail Here4 itself as the CAN-to-PWM converter (via its second breakout cable) plus a custom tail PCB for regulation. That plan is dropped. Here4 doing both secondary GPS and tail conversion would let a single Here4 failure remove tail control-surface outputs, not just GPS, and Here4/RTK failures have been seen on Quiver, so a dedicated adapter avoids that single point of failure. Decoupling the servo signal path from the GPS also removes the dependency on the Here4 breakout, and an off-the-shelf adapter removes the custom PCB spin. The former Here4 breakout PWM bench test is no longer needed.

#### CAN-to-PWM Adapter — Matek CAN-L4-PWM

| Item | Detail |
|---|---|
| Function | DroneCAN node. Receives RuddervatorLeft/Right (SERVO_FUNCTION 79/80) over CAN, outputs PWM |
| CAN | 4-pin JST-GH, daisy-chained from CAN1. Node ID 12, 120Ω term if at bus end |
| Logic power | 4.5–5.5V from the tail PM12S-3 (not the CAN 5V rail), ~30mA |
| PWM outputs | 9× PWM (8× DShot-capable), use 2 for ruddervator L/R at 5V logic. AP_Periph firmware (STM32L431), 3.5g |
| Servo rail | V+ rail fed from the tail PM12S-3 Vx output (5.25V PT1, 8V final). Adapter logic and Here4 stay on the PM12S-3 fixed 5V, deliberately separate so the later 8V switch cannot reach the logic or the GPS |
| ArduPilot | Map the CAN node's PWM channels to SERVO_FUNCTION 79/80, confirm in WP-E04 |

Sourcing (June 11): in stock in Turkey at Aykut Havacılık, 2,060 TL. The Matek CAN-L431 from the same shop (2,500 TL, 5× PWM, 3× UART, I2C, dual parallel GH-4P CAN connectors) is an acceptable alternate running the same AP_Periph family on the same STM32L431. The L4-PWM stays primary: more PWM headroom, and it sits at the bus end where its single CAN connector plus the 120Ω terminator fits the topology. The L431's pass-through CAN connectors would matter only if the adapter had to sit mid-bus.

#### Tail Servo Power — Matek PM12S-3

| Item | Detail |
|---|---|
| Input | HV bus spur (36–50V), an 18 AWG branch broken out at the tail-boom ESC-lead extension joint (not a dedicated FCHUB run), via the tail harness |
| Outputs | Vx at 5.25V (PT1) / 8V (final), 15A cont, for the 2 ruddervator servos. Fixed 5.2V/4A for the adapter logic and the Here4. Fixed 12V unused at the tail |
| Type | Matek PM12S-3 (9–55V input with TVS, 53g). The 8V switch for the final 6–8.4V servos is a pad setting, no hardware change |
| Wiring | Vx lands on the adapter's servo V+ rail, servos plug into the adapter PWM headers. The fixed 5V feeds the adapter logic pad and the Here4 on a separate feed |
| Kill scope | On the HV spur, de-energized only by the HV kill, not the LV kill (see §4.3) |

Common ground: all PM12S-3 rails share ground, so the Vx servo rail, the adapter logic, the Here4, and the PWM signals reference one ground.

**Kill behavior:** the tail PM12S-3 feeds the servos, the adapter, and the Here4 from the HV spur, so the entire tail control path is on the HV kill only and is unaffected by the LV kill. There is no longer a dependency on the CAN 5V rail for tail control.

**Fallback for first hover:** for the first VTOL-only bench/hover test, the adapter and servos can be powered from any bench 5V source before the tail harness and PM12S-3 are finalized.

#### Wing Flaperon Drive — 6X IO PWM (wing CAN node considered and dropped)

A second CAN-L4-PWM at the wing center was considered on June 11 and dropped the same day with the Pixhawk 6X selection: the 6X IO PWM outputs 1–2 drive the flaperons directly, saving a node, its wiring, and shelf space. The flaperon PWM runs (~1–1.5 m to the wing roots) follow the §4.6 EMI rules. Servo power is unchanged: the forward PM12S-3 Vx feeds the flaperon V+ through the wing connectors. Only the tail keeps a CAN-PWM adapter (node 12), so `CAN_D1_UC_SRV_BM = 0x000C` (servo outputs 3–4).

---

### 4.11 Servo System

**PT1 servo voltage (June 2026):** PT1 uses 5V test servos for both wing and tail. The first prototype only needs hover and control-surface motion, not loaded cruise, so 5V servos off the existing 5V rails simplify the wiring and drop the HV servo rails. This is a test choice, not the final loaded selection. The final design is expected to use 6–8.4V HV servos sized to Zeynep's control-surface hinge-moment study (WP-E06), which switches both PM12S-3 Vx rails (tail and forward) to 8V with no hardware change. The torque and current targets below are the final-design requirements, pending Alperen's main-wing and tail numbers.

#### Ruddervator Servos (×2)

V-tail airfoil: NACA 0015, 280mm chord, control surface 30–35% chord (~84–98 mm deep surface). Each surface is one ruddervator providing mixed pitch + yaw authority (ArduPilot SERVO_FUNCTION 79/80).

| Parameter | Requirement |
|---|---|
| Type | PT1: 5V digital. Final: digital HV (6–8.4V) |
| Torque target | ≥ 8 kg·cm at 7.4V |
| Speed | ≤ 0.12 s / 60° at 7.4V |
| Weight | ≤ 30g each (tail weight is CG-critical) |
| Connector | Standard 3-pin JR/Futaba (2.54mm pitch) |
| Control | Via tail CAN-PWM adapter PWM output (see §4.10) |

Torque target is preliminary — Zeynep's control surface sizing study output (WP-E06 input) will set the final requirement. Verify hinge moment calculation before locking servo selection.

#### Flaperon Servos (×2)

Main wing: Clark Z, 20–25% chord flaperons.

| Parameter | Requirement |
|---|---|
| Type | PT1: 5V digital. Final: digital HV (6–8.4V) |
| Torque target | ≥ 6 kg·cm at 7.4V |
| Speed | ≤ 0.15 s / 60° at 7.4V |
| Weight | ≤ 30g each |
| Connector | Standard 3-pin + wing panel connector (see §4.6) |
| Control | Via Pixhawk 6X IO PWM outputs 1 and 2 (§4.6 EMI rules on the wing runs) |

These are final-design HV candidates: Savöx SH-0255MG (~20g, 9 kg·cm at 7.4V), KST DS135MG (~27g, 13 kg·cm at 7.4V), MKS HV6130 (~20g, 10 kg·cm at 7.4V). Zeynep also flagged a Turkey-available servo with a 6–8.4V range for the final design. For PT1, use inexpensive 5V digital servos (hover and control-surface motion only). Selection in WP-E06.

---

### 4.12 Telemetry and Communications

#### RC Link

| Parameter | Value |
|---|---|
| Range requirement (Phase 1) | VLOS, ≤ 5 km |
| Protocol | CRSF (preferred, bidirectional) or SBUS |
| FC interface | UART5 (SERIAL5) |
| Receiver location | Avionics shelf or fuselage side, antenna external |
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

The pack is a dumb LiPo, so monitoring is sensor-based rather than from a BMS. Vehicle current comes from the FCHUB-12S 440A analog sensor, pack voltage from a FC voltage sense or the power module (analog on the PM02 V3, I2C on the PM02D HV), and SoC from ArduPilot coulomb counting. There is no per-cell data.

| Parameter | Method |
|---|---|
| Pack current | FCHUB-12S 440A analog sensor (3.3V ADC, scale 75) → FC ADC pin |
| Pack voltage | FCHUB voltage sense → FC ADC, or the power module (analog on 6C, I2C on 6X) |
| Cell voltages | Not available (no BMS). Use a plug-in cell checker/alarm for ground checks |
| SoC / consumed | ArduPilot coulomb counting from the FCHUB current sensor |
| Temperature | Pack not instrumented (no BMS). Optional external thermistor if needed |
| Low voltage warning | BATT_LOW_VOLT = 42V (3.5V/cell), GCS alert + buzzer |
| Critical voltage | BATT_CRT_VOLT = 39.6V (3.3V/cell), RTL trigger |
| Low capacity | BATT_LOW_MAH = 2000 mAh (reserve threshold) |

ArduPilot parameter: `BATT_MONITOR = 4` (Analog Voltage and Current). Set `BATT_CURR_PIN` / `BATT_VOLT_PIN` to the Pixhawk ADC pins wired to the FCHUB sense outputs (the 6X keeps ADC pins on its AD&IO port) and `BATT_CURR_SCALE` per the FCHUB (3.3V ADC, scale 75). Voltage sense is confirmed on the FCHUB-12S V2: the V pad set to VBat feeds through the built-in 1K:20K divider, so 50.4V full charge reads ~2.4V at the ADC and `BATT_VOLT_MULT` ≈ 21 (calibrate against a multimeter). The power module stays available as BATT2.

**Power module as secondary monitor:** the module gives a redundant pack-voltage reading as BATT2 (analog pins on the PM02 V3, I2C on the PM02D HV). Its current reflects only the FC branch, not vehicle current, because the module feeds a dead-end branch off the FCHUB pad tap. Vehicle current and consumed-mAh come from the FCHUB sensor (primary).

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

Updated June 11 for Alperen's fuselage conceptual v5: longer thinner nose, battery moved forward for CG, avionics on a large shelf above the battery with top access, pusher between the fuselage and the boom-carried tail. Detailed layout remains WP-E07, pending shelf and bay dimensions.

```
Nose (slimmed, v5):
  - No avionics requirement. Pitot mast and antenna placement only if needed

Forward battery bay:
  - ProFuse 12S pack (186×76×160mm, 3,822g): strapped, fore-aft position sets CG
  - XT150 pack connector: reachable for ground safing
  - HV kill (OQ-05): adjacent, actuation per the chosen implementation

Avionics shelf (above battery, top access):
  - Pixhawk 6C/6X: vibration-isolated, X-axis forward
  - Power module (PM02 V3 or PM02D HV): inline on the pad-tap leg
  - Primary GPS: on the shelf cover, unobstructed sky view
  - RC receiver + telemetry radio: shelf edges, antennas external
  - LV kill switch: reachable through the shelf hatch
  - Shelf footprint reserves rectangular space for a future integration PCB (per Alperen)

CG bay / boom roots:
  - FCHUB-12S PDB: near the boom-root convergence, short ESC lead runs
  - CAN splitter/hub: adjacent to the FCHUB

Booms (×4 motor arms, 2 carrying the tail):
  - ESCs: externally mounted near motors (within 150mm of motor centerline)
  - Phase wire exits sealed where they enter the boom
  - Tail harness (CAN trunk + 18 AWG HV spur off the ESC-lead extension joint): along one tail boom
  - Inside-tube routing is attractive but waits on the structural check (June 5)

Tail group (boom-mounted):
  - Here4 (secondary GPS): top surface, unobstructed sky view
  - CAN-PWM adapter + PM12S-3: bulkhead-mounted, connectors accessible
  - Ruddervator servos: local, short leads

Wings (detachable panels):
  - Flaperon servo: inline with control horn, servo output axis on hinge line
  - Wing panel connector: at wing root, mates to fuselage connector on insertion
```

---

### 4.16 Phase 1 BOM and Cost Estimate

Procurement view with priorities, candidate vendors, and lead flags: `docs/pt1-parts-list.md` (SPH-E-003). The tables below stay the cost baseline.

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
| Pixhawk 6X Standard Set v2A (FC + baseboard + PM02D HV + M9N GPS + cables, selected June 11) | 1 | $321 | $321 |
| ProFuse Super Nano 12S 16Ah 60C (XT150, robotsepeti TR, 41,322 TL) | 1 | ~$1,030 | ~$1,030 |
| Matek FCHUB-12S V2 PDB (HV distribution + 5V/12V BECs), incl. 1 spare | 2 | $30 | $60 |
| XT60 pigtail + 16 AWG for the FCHUB pad tap (Y-splitter dropped) | 1 lot | $10 | $10 |
| Primary GPS, nose (DroneCAN, model TBD, OQ-12) | 1 | TBD | TBD |
| Here4 GPS, tail (secondary GPS, GPS only) | 1 | $300 | $300 |
| Matek PM12S-3 ×2 (tail + forward servo rails, Aykut Havacılık TR, 6,900 TL each) | 2 | $170 | $340 |
| HV kill hardware (D2 outcome: EV200 ~€149, MOSFET switch ~$80, or software-only $0) | 1 | TBD | TBD |
| LV kill switch | 1 | $15 | $15 |
| XT150 male connectors (pack mate) + optional AS150U adapter pigtail | 1 lot | $15 | $15 |
| XT90S connector pairs (ESC boom junctions + spares) | 8 | $6 | $48 |
| 5mm bullet connectors (×12 sets) | 12 | $2 | $24 |
| 8 AWG silicone wire (2 m, battery harness) | 1 lot | $15 | $15 |
| 8 AWG silicone wire (10 m) | 1 lot | $20 | $20 |
| 12 AWG silicone wire (extension leads) | 1 lot | $12 | $12 |
| LV wiring assortment (18–26 AWG) | 1 lot | $30 | $30 |
| CAN bus twisted-pair wire (5 m) | 1 lot | $10 | $10 |
| JST-GH crimp kit (CAN, UART, assorted) | 1 | $25 | $25 |
| CAN splitter/hub (DroneCAN, JST-GH, model TBD) | 1 | $15 | $15 |
| 5V digital ruddervator servos (PT1) | 2 | $25 | $50 |
| 5V digital flaperon servos (PT1) | 2 | $25 | $50 |
| Wing panel servo connectors (Molex Mini-Fit Jr) | 2 sets | $10 | $20 |
| RC receiver (TBD — pending RC system decision) | 1 | TBD | TBD |
| Telemetry radio (TBD — pending SIYI HM30 decision) | 1 | TBD | TBD |
| Ainstein US-D1 radar altimeter (selected June 11, believed on hand in TR, confirm) | 1 | $0–280 | ~$0 |
| Pitot tube transducer (if needed beyond existing) | 1 | $30 | $30 |
| WS2812B LED strip + wiring | 1 lot | $20 | $20 |
| Matek CAN-L4-PWM adapter (tail), incl. 1 spare (Aykut Havacılık TR, 2,060 TL each) | 2 | $50 | $100 |
| (tail and flaperon servo rails covered by the PM12S-3 line above) | — | — | — |
| Heat-shrink, zip ties, mounting hardware | 1 lot | $40 | $40 |
| **Subtotal (to procure, est.)** | | | **~$2,550 + TBD items** |
| **Phase 1 Grand Total (est., excl. labor)** | | | **~$5,185 + TBD items** |

Estimates only. Verify against actual quotes before procurement. Some items may be on hand from Quiver or other projects.

---

## 5. Phase 2: Hybrid Integration

Target: July–October 2026. Objectives: IC engine integration, transition flights, vibration characterization and notch filter tuning. Engine model for current layout work is a DLE55cc-class gasoline engine (June 11). The generator path moved out of Phase 2: PT2 may simply fly with a gasoline engine, and starter-generator work is now a PT3 / end-of-year sub-project candidate. Three starter-generator vendors responded to Alperen (Czech-linked, German, US), with a US vendor call pending. The §5.8 VESC architecture below predates this and reads as background for that research.

> **⚠️ STATUS: This section is preliminary. Significant research is still required before any of the Phase 2 content below should be treated as a design baseline.**
>
> The IC engine, fuel system, starter, generator, ignition, RPM sensor, throttle servo, and vibration isolation subsections in §5 are working notes based on general 2-stroke RC engine practice and a candidate engine class (35–55cc Stinger / DLE). The actual engine has not been selected or purchased, no bench testing has occurred, and most parameters (fuel consumption rate, electrical load profile, vibration spectrum, alternator output, starter inrush current) will only be known after empirical measurement on the chosen unit. Treat all numbers, parts, and architectures in §5 as candidate values to be validated, not commitments. Specific items requiring further research before commitment:
>
> - Engine selection (Stinger vs DLE, 35cc vs 55cc) — buy both and test empirically per current plan
> - Starter motor selection and current draw (Pilot RC: ~30 days lead after order, and the vendor says it cannot double as a generator)
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

The Pilot RC DLE35RA electric auto-start unit (420g) is a gear-driven brushless motor that engages the engine crankshaft for starting. Vendor status (May 2026): not in stock, can be prepared roughly 30 days after order. The vendor stated the unit cannot be used as a generator (reason not given), so a generator path would need a separate VESC and brushless motor (§5.8), not this starter.

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

### WP-E03 — Tail CAN-PWM Adapter and Servo Power Integration

| Field | Value |
|---|---|
| Owner | Erick |
| Phase | 1 |
| Inputs | Matek CAN-L4-PWM adapter, PM12S-3 (Vx servo rail + fixed 5V logic feed), tail bay dimensions from Alperen, servo connector type from WP-E06 |
| Outputs | Mounted and wired tail electronics: HV spur → PM12S-3 → adapter servo rail (Vx) + logic/Here4 feed (5V), CAN daisy-chain, servo harness |
| Dependencies | WP-E04 (adapter PWM validation); WP-E06 for servo connector type; Alperen's tail bay dimensions |
| Deliverables | Tail wiring diagram, harness, bench-verified servo response, mounting solution. No custom PCB (off-the-shelf adapter + PM12S-3). |
| Note | Supersedes the earlier custom tail PCB. Revisit a custom integration board only if mounting two loose modules proves impractical. |

### WP-E04 — CAN-PWM Adapter Validation

| Field | Value |
|---|---|
| Owner | Erick (executor), Zeynep (ArduPilot config) |
| Phase | 1 |
| Inputs | Matek CAN-L4-PWM, ArduPilot-configured FC (existing Quiver hardware acceptable), oscilloscope or servo response check |
| Outputs | Pass/fail: adapter enumerates as a DroneCAN node, accepts SERVO_FUNCTION 79/80 mapping, outputs correct PWM timing at 5V logic, servo responds |
| Dependencies | Adapter in hand; CAN bus config from §4.9 |
| Deliverables | Short test report (DroneCAN node ID, parameter mapping, scope/servo result), posted to Spearhead forum |
| Note | Replaces the former Here4 breakout PWM test. Confirm ArduPilot 4.x DroneCAN servo output support for the chosen adapter firmware. |

### WP-E05 — Avionics and Sensor Integration

| Field | Value |
|---|---|
| Owner | Erick (integration), Zeynep (ArduPilot configuration) |
| Phase | 1 |
| Inputs | Pixhawk 6C pinout, GPS CAN pinouts (nose module + tail Here4), altimeter selection (Ainstein vs Benewake), pitot tube model from storage, CAN bus topology from §4.9 |
| Outputs | Full avionics wiring diagram, ArduPilot parameter file (.param) covering all sensors, validated bench-tested avionics stack |
| Dependencies | WP-E04 (CAN-PWM adapter confirmed); altimeter selection is part of this WP |
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
| Inputs | Component dimensions (Pixhawk 6C/6X, BECs, power module, RC receiver, telemetry radio), avionics shelf dimensions from Alperen (fuselage v5), forward battery bay dimensions, CG target from Alperen |
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
| OQ-01 | Battery: pack format (one 12S vs 2× 6S in series), brand, and forward-bay CG position locked | Erick / Alperen | High | All Phase 1 work |
| OQ-02 | CAN-PWM adapter (Matek CAN-L4-PWM): enumerates as a DroneCAN servo node and outputs correct PWM under ArduPilot 4.x? | Erick / Zeynep | High | WP-E04 |
| OQ-03 | ~~Altimeter: Ainstein US-D1 vs Benewake TF03-180?~~ Closed June 11: US-D1 selected (radar robustness in prop wash dust). Quote pending, longest lead item | Erick | — | — |
| OQ-04 | RC system: ELRS vs Crossfire? (depends on team transmitter) | Alperen / Erick | Medium | WP-E08 |
| OQ-05 | HV kill: must be RC-switchable (June 11 criterion). EV200 contactor + relay driver vs MOSFET anti-spark switch vs software E-stop + recorded deviation (§4.3) | Erick / Alperen | High | WP-E02, PT1 order |
| OQ-06 | Pitot tube (from storage): model and interface confirmed? | Alperen | Medium | WP-E05 |
| OQ-07 | Shared or separate battery for VTOL vs IC ignition in Phase 2? | Erick | Low | Phase 2 design |
| OQ-08 | Electric starter generator mode: what is the mechanical limitation? | Erick | Low | WP-E11 |
| OQ-09 | Wing panel connector: Molex Mini-Fit Jr confirmed, or alternative? | Erick | Medium | WP-E06 |
| OQ-10 | Payload quick-release connector: align with Arrow platform standard? | Erick / Alperen | Low | Phase 4 design |
| OQ-11 | ~~Here4 breakout board: does Thomas's kit include one?~~ Closed, superseded by the standalone adapter. Here4 no longer used as converter. | Thomas | — | — |
| OQ-12 | Primary GPS: second Here4 (matched RTK pair, enables GPS yaw) or F9P-class/salvage with on-hand antennas (compass yaw)? Decide before the PT1 order ships (parts list D1) | Team | High | WP-E05, PT1 order |

---

## 10. Revision History

| Rev | Date | Author | Changes |
|---|---|---|---|
| 0.16 | June 11, 2026 | errrks | Reconciled with the June 5 and June 11 calls. FC: Pixhawk 6C or 6X by sourcing (Zeynep approved either). Power module follows the FC and the pairing error is fixed: PM02D is I2C, 5X/6X only, so the 6C path uses the analog PM02 V3 and the 6X set includes the PM02D HV. Battery Y-splitter dropped: the FC power module taps the FCHUB battery input pads (Y stays as fallback). Component zones reworked for fuselage v5: forward battery bay, top-access avionics shelf, slimmed nose, boom-routed tail harness. HV kill: added the RC-switchable criterion (June 11) and a new candidate table (EV200 + relay driver, MOSFET anti-spark switch, software E-stop + recorded deviation), superseding the PP75 option (undersized at 75A class vs 98–124A hover). Altimeter selected: US-D1, specs corrected to 110g / 100 Hz / UART+CAN per Ainstein. Regulators standardized on Matek BEC12S-PRO (5.2/8/12V selectable). Added the PT1 harness connection table (§4.6, H1–H21), renamed §4.6, and split the procurement view into `docs/pt1-parts-list.md` (SPH-E-003). §5: DLE55-class layout model, starter-generator moved to PT3 research. Late June 11 additions: battery selected (ProFuse Super Nano 12S 16Ah 60C, XT150 socket, 3,822g, 186×76×160mm, robotsepeti TR), battery interface moved from AS150U to XT150 with the anti-spark burden shifted to the inline kill (§4.3 connector interaction note), FCHUB-12S confirmed as V2 (built-in 1K:20K VBat divider closes the §4.13 voltage sense question, 12VSW PINIO pad can drive the EV200 coil directly), tail CAN boards sourced in Turkey at Aykut Havacılık with the CAN-L431 documented as alternate, US-D1 believed on hand in Turkey pending confirmation. Regulators changed to 2× Matek PM12S-3 (Vx 15A selectable 5.25/6/8V + fixed 5V + fixed 12V per unit, Aykut Havacılık), superseding the BEC12S-PRO plan: tail unit feeds the ruddervator Vx rail plus adapter/Here4 5V, forward unit feeds the flaperon Vx rail, and the final 8V servo change becomes a pad setting. A wing CAN-PWM node (13) was added and dropped the same evening with the FC decision. Final state: **Pixhawk 6X selected (D4 closed)**, chosen because its IO PWM outputs 1–2 drive the wing flaperons directly, saving the second adapter and shelf space. CAN-L4-PWM order stays at 2 (tail node 12 + spare), `CAN_D1_UC_SRV_BM = 0x000C`, and the §4.6 EMI rules re-apply to the ~1–1.5 m wing PWM runs. 6X sourcing: Aykut lists the set at 38,220 TL (~3× Holybro direct) and is out of stock, so Holybro EU at $321 is the order path. The set's M9N GPS is kept as a bench/backup unit. PM12S-3 TR markup vs Matek direct flagged in the parts list. Battery harness corrected from 4 AWG to 8 AWG: 4 AWG does not fit XT150 solder cups and the connector is the thermal limit on the short run. |
| 0.15 | June 2026 | errrks | Reconciled with the June 1 call (PT1 decisions). Servos: PT1 uses 5V test servos for wing and tail, dropping both 7.8V BECs. Flaperons run on the peripheral 5V rail, the tail runs on a local HV→5V BEC that also powers the adapter and Here4, which resolves the tail-5V run and puts the whole tail control path on the HV kill only. Final design keeps 6–8.4V HV servos. Battery splitter reverted from XT90 to AS150U/XT150 (full current, 4AWG). Added a §4.6 routing and EMI subsection. Kill switch: PT1 implementation left TBD (contactor/SSR vs manual disconnect), §9.3 requirement noted. Updated §2.1/§2.3, §4.1, §4.3, §4.4, §4.6, §4.10, §4.11, and the BOM. |
| 0.14 | June 2026 | errrks | Battery changed from a DroneCAN smart pack to a dumb 12S LiPo at 16Ah (Alperen's May mission-energy rerun, saves ~1kg). Monitoring moves to the FCHUB-12S 440A analog current sensor plus a voltage sense (BATT_MONITOR = 4), with the PM02D I2C as a redundant voltage reading. No per-cell or CAN battery data. Reworked §4.1, §4.2, §4.13, §2.3, and the BOM. Reconciled against the May 2026 meeting notes. |
| 0.13 | June 2026 | errrks | Flight controller powered by a Holybro PM02D (XT60 input, 5.2V/3A + I2C) on its own leg of an XT90 splitter at the battery, isolated from the PDB. Single FC supply: POWER1 = PM02D, POWER2 unused (no hot-standby for V1). Restructured §4.4 5V rails into FC power, peripheral 5V, and LED 5V. XT90 accepted for the high-current path (short hover duty, monitor connector temps), high-current leads drop to ~10AWG to fit XT90, the PM02D leg is a smaller-gauge XT60 branch. PM02D I2C is an optional BATT2 voltage monitor and its current reads only the FC branch. Updated §2.1/§2.3, §4.3, §4.4, §4.6, §4.7, §4.13, and the BOM. |
| 0.12 | June 2026 | errrks | Power distribution: selected the off-the-shelf Matek FCHUB-12S PDB for V1 (no custom PDB). Corrected its specs (4×70A continuous and 4×110A burst per ESC pad, 8–60V input, on-board 5V/5A + 12V/4A + 3.3V BECs, 440A current sensor left unused, solder pads not screw terminals). Only the power pads are used since the ESCs run DroneCAN. CAN1 to the ESCs now fans out through a passive CAN splitter/hub instead of a custom board. Updated §2.1/§2.2/§2.3, §4.2, §4.4, §4.9, and the BOM. |
| 0.11 | June 2026 | errrks | Tail architecture change: the tail Here4 is now GPS only and a standalone DroneCAN-to-PWM adapter (Matek CAN-L4-PWM, node 12) drives the ruddervators. Custom tail PCB dropped. Ruddervator 7.8V now comes from an off-the-shelf HV→7.8V BEC on the HV spur. Confirmed two redundant GPS (AIP-006 §9.4): primary nose module (DroneCAN, TBD) plus secondary tail Here4. Rewrote §4.10 and updated §2.1/§2.2/§2.3, §4.3/§4.4/§4.6/§4.7/§4.8/§4.9, the BOM, WP-E03/E04/E05, and OQ-02/OQ-11. Added OQ-12. Note: a mixed GPS pair does not support moving baseline yaw. |
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
