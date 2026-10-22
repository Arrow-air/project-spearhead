---
title: "Spearhead PT2 PCB Scope and Reuse Map"
sidebar_label: "PT2 PCB Scope and Reuse Map"
sidebar_position: 10
---

# Spearhead PT2 PCB Scope and Reuse Map

**Doc ID:** SPH-E-008 (draft, renumbered October 16, 2026, was SPH-E-005)
**Revision:** 0.4, September 28, 2026
**Author:** errrks
**Status:** Draft for the September 18 team call. Nothing here is a decision until the call closes it.

Companion to SPH-E-001 (electrical master, Rev 0.19), SPH-E-002 (PT1 decisions), SPH-E-003 (PT1 parts list).

---

## 1. Purpose

Define the PT2 PCB set for Spearhead by reusing the Quiver PT3 and Caribou boards, and list what has to be decided before layout starts. Design target is the November 2026 US build. Design window is 2 to 3 weeks, then roughly one month for fab, parts and shipping (lead time to be confirmed with Julius).

Source of the scope: the September 15 call (Alperen, Erick, Thomas, Zeynep, KBM). Source of the reuse inventory: the KiCad files listed in section 3, read on September 18.

---

## 2. What PT1 taught us

PT1 flew on the bare minimum. Verified from build photos, not from any document:

| SPH-E-001 section | Planned for PT1 | Flown on PT1 |
|---|---|---|
| 4.2, 4.13 FCHUB with two-port sensing | Vehicle V and I into POWER2 as BATT1 | PDB only, no sensing |
| 4.4 Regulation | Wing buck at 8.4 V, tail PM12S-3 | None |
| 4.5, 4.9 MAD DroneCAN ESCs on booms | CAN nodes 1 to 4 plus tail 11 and 12 | KDE PWM ESCs in the fuselage, nothing on CAN |
| 4.10, 4.11 Tail node, Here4, four Kingmax servos | Full control surface path | Nothing wired, wings off |
| 4.14 LEDs, 4.8 altimeter and pitot | Included or optional | Absent |

Thomas's notes from the build and first flight (September 6) are the second as-built source and agree with the photos:

- 12S battery into the PDB (he wrote "XT90 in?", the Aug 27 decision was QS8 on the pack side, XT90 female on the PDB live side).
- PDB to ESCs 1 to 4. KDE ESCs in the cockpit.
- PDB to a 5 V converter for FC and servo power, with opto isolation.
- PM02 power module for the Pixhawk. Pixhawk powers all current LV parts.
- Here4 GPS, telemetry radio and RC receiver all in the cockpit. Nothing in the tail.

Consequences for PT2:

- The CAN architecture, tail architecture, both servo rails, and battery monitoring have never flown. Zeynep's June bench test is the only validation of the tail node.
- The PT1 PDB has no spare HV output. PT1.5 transition uses a separate pusher-only battery.
- The bay was six loose boards and cabling. Cable volume, not component volume, drove the packing problem.

---

## 3. Source boards (verified from KiCad files, September 18)

| Board | Date | Outline | Cu layers | Footprints | What it is |
|---|---|---|---|---|---|
| Quiver Front (battery) PCB V1.3 | Jan 5, 2026 | 247.5 x 167.5 mm | 4 | 87 | HV path from the Tattu pack. Custom Molex 46437 battery mate, bolt-in fuse, five parallel N-MOSFETs (S400N08TL) as the main switch, precharge with bypass, two CPC1019N SSRs, TSR 48 V to 12 V regulator, TSYS01 temperature sensor, CAN and I2C to the main board, status LEDs, Amphenol AMT0450003 screw terminals |
| Quiver Main PCB V1.2 | Jan 5, 2026 | 356.0 x 253.7 mm | 6 | 164 | Avionics hub. Pix32 V6 carrier interface (2x AXK6S mezzanine), RPi 5 header, two Gigablox Nano Ethernet switch modules, PPP-to-Ethernet adapter ports, MCP2515 and TJA1049 CAN for the Pi, REC20K 5 V and REC30K 12 V isolated converters, three CPC1019N SSRs, SIRA99DP and SIR570DP MOSFET switches, CAN2 termination switch, TSYS01 temperature sensor, buzzer, four XT60PW, one XT30PW, M9N GPS footprint, F9P ports, USB-C, two M6 power taps, Littelfuse 0154 fuses and a PTC |
| Quiver FC PCB V1 | Sep 29, 2025 | 200.1 x 134.0 mm | 4 | 52 | Pix32 V6 carrier. AXK5S mezzanine connectors J1 to J3, breakouts for PWM, GPIO, GPS, I2C, CAN, SBUS and power per the Holybro schematic |
| Caribou CMAIN | Jun 30, 2026 | 331.3 x 260.0 mm | 4 | 248 | Quiver Main merged with the Feather PDB for a six-motor frame. Same Pi, Ethernet, CAN and mezzanine blocks. Adds three LTC4370 ideal-diode 12 V combiners, three SY8205 bucks, nine CPC1019N SSRs, JST JWPF waterproof connectors throughout, twelve status LEDs. Has a documented net-class and routing-rule set |
| Caribou CBC | Jun 30, 2026 | 265.5 x 290.5 mm | 6 | 301 | Per-arm 18S battery board. Precharge, hardware kill latch, main enable on TOLL-8 MOSFETs, isolated 12 V at 10 A push-pull, 5 V and 3.3 V logic, dual CAN (battery protocol and DroneCAN bridge, ADM3053 isolated transceiver), ESP32-S3 WiFi and BLE diagnostics, 250 A bolt-in fuse |
| Caribou C90D | Jun 9, 2026 | 150 x 30 mm | 6 | 10 | 90 degree battery adapter. Not relevant to Spearhead |

File locations (latest known copies):

- Quiver Front V1.3: `~/Desktop/Arrow/2026/1-5-26/Front_PCB_V1_3/Front_PCB/Front_PCB.kicad_pcb`
- Quiver Main V1.2: `~/Desktop/Arrow/2026/1-5-26/Quiver_PT3_Main_PCB_V1_2/Quiver_PT3_Main_PCB/Quiver_PT3_Main_PCB.kicad_pcb`
- Quiver FC PCB V1: `project-quiver/task-grant-bounty/pt3/electronics/0008-FC-PCB/KiCAD-Files/Quiver_PT3_FC_PCB.kicad_pcb` (identical to the 9-29-25 Desktop copy)
- Caribou boards: `project-caribou/engineering/electronics/pcbs/{CMAIN_PCB,CBC_PCB,C90D_PCB}/kicad/`
- Quiver information notes: `project-quiver/task-grant-bounty/pt3/electronics/000{1,2,7,8}-*/`

Notes on the files:

- Quiver Main V1.2 has 6 copper layers in the file. The V1 information note says 4. Confirm which stackup was fabricated before copying it.
- Quiver boards use Samacsys, Snapeda and EasyEDA library parts that are not stored in the project. Caribou boards keep libraries local under `kicad/libs/`. Spearhead should follow the Caribou convention from day one.
- Caribou CMAIN is a KiCad 10 project. Quiver boards are KiCad 9 format. Open everything in KiCad 10 and let it upgrade.

---

## 4. Proposed Spearhead board set

Quiver is three boards: battery, main, and FC carrier. Caribou is a per-arm battery board plus a main board. Spearhead maps as follows.

| Spearhead board | Role | Primary source | Location |
|---|---|---|---|
| **SBAT, battery board** | HV entry, fuse, hardware kill, precharge, current and voltage sense, pack and board temperature, HV outputs to ESCs and pusher | Quiver Front PCB V1.3 | Forward battery bay, next to the pack |
| **SMAIN, cockpit board** | Pix32 V6 carrier, regulation, CAN fan-out, Pi and Ethernet, peripheral ports, servo rails, LV kill | Quiver Main V1.2 plus Quiver FC PCB, Caribou CMAIN for the SY8205 bucks and JWPF connectors | Avionics shelf |
| **STAIL, tail board** | HV in, 8 V and 5 V rails, CAN-L4-PWM hat, Here4 connector, two servo headers | New, small. Regulator block borrowed from Caribou CBC 5 V buck or Quiver TSR | Tail bay |
| ~~SPUSH, pusher hat~~ | **Dropped September 25, 2026.** The electric pusher runs from a secondary Longshot battery in the tank bay, wired directly to the pusher ESC. SMAIN carries a permanent port set for both variants instead: throttle PWM, engine temperature I2C, ignition power output, pusher ESC CAN. See `engineering/electronics/pcbs/SMAIN_PCB/docs/pinout.md` A6 and section 5.4 below | | |

Why a separate SBAT rather than folding HV into SMAIN:

- It restores the hardware HV kill that AIP-006 section 9.3 requires and PT1 deviated on (SPH-E-001 section 4.3, SPH-E-002 deviation register). The Quiver Front PCB already implements exactly that: MOSFET main switch, precharge, SSR control from the FC, and a pushbutton precharge start.
- It keeps 320 A peak copper off the avionics board, so SMAIN can be thin and dense.
- The battery bay and the avionics shelf are different compartments in Alperen's fuselage v5 layout. A single board would span both.

### 4.1 PT2 system architecture

Power topology. Every line that crosses a board boundary is a harness.

```
                 Longshot 14S3P ~19.5 Ah pack, AS150U (4 signal pins), ~4 kg, 58.8 V full
                                   |
                                   | 8 AWG
                 +-----------------v------------------+
                 |  SBAT  (forward battery bay)        |
                 |  fuse -> precharge -> MOSFET kill   |  <- SSR enable + precharge end from SMAIN (2 signals)
                 |  shunt or hall current sense        |  -> V and I analog to SMAIN (Pix32 power sense)
                 |  VBat divider, TSYS01 temp (I2C)    |  -> I2C to SMAIN
                 |  status LEDs, precharge pushbutton  |
                 +--+------+------+------+------+---+--+
                    |      |      |      |      |   |
                 ESC1   ESC2   ESC3   ESC4   pusher  HV spur, 18 AWG
                 12 AWG each, 80 A class     HV out  |
                 (ESCs mounted next to SBAT) (to SPUSH)   |
                                                       v
      +--------------------------------------+   +-----v---------------------+
      |  SMAIN  (avionics shelf)              |   |  STAIL  (tail bay)         |
      |  HV in (16 AWG from SBAT) ->          |   |  HV in -> 8 V Vx 15 A      |
      |    5 V and 12 V rails (fused, LEDs)   |   |          -> 5 V logic      |
      |    8.4 V wing servo rail, 15 A        |   |  CAN in/out, 120 R jumper  |
      |    LV kill input (FC 5 V + periph 5 V)|   |  CAN-L4-PWM hat            |
      |  FC PCB stack on J3, J4 sockets       |   |  Here4 connector (cabled)  |
      |  Pi 5 + Ethernet switch + PPP bridge  |   |  2x ruddervator JR headers |
      |  CAN1 fan-out, CAN2 reserved          |   +----------------------------+
      |  GPS1 serial (F9P), TELEM1, RC IN,    |
      |  UART for US-D1, I2C for pitot/temps, |         +----------------------+
      |  LED strip feed + NeoPixel GPIO       |         |  SPUSH  (aft, hat)   |
      |  2x wing connectors (PWM, 8.4 V, GND) |         |  electric: HV out,   |
      |  payload 12VSW, SIYI ports, USB-C     |         |   ESC signal + telem |
      +---------------------------------------+         |  gasoline: CDI supply|
                                                        |   4.8-8.4 V isolated,|
                                                        |   manual switch +    |
                                                        |   opto kill, tach in,|
                                                        |   throttle/choke PWM |
                                                        +----------------------+
```

Communications topology.

```
Pix32 V6 (on the FC PCB, stacked on SMAIN)
|
+-- CAN1 --[120R on SMAIN]-- ESC1..4 stubs (DroneCAN, MAD) -- trunk to STAIL -- Here4 -- CAN-L4-PWM --[120R jumper]
|                            (each ESC connector also carries a PWM pin for the KDE fallback)
+-- CAN2 --[term switch]-- reserved: Remote ID node, LED node, payload
+-- SPI/UART/CAN to Pi 5 (MCP2515 + TJA1049 block), Pi Ethernet -> Gigablox switch -> PPP bridge, SIYI, payload
+-- GPS1 serial -- F9P (shelf cover)
+-- TELEM1 -- SiK V3 now, RFD900x class later
+-- RC IN -- R9DS SBUS
+-- UART -- US-D1 altimeter
+-- I2C -- pitot transducer, TSYS01 on SMAIN, TSYS01 on SBAT
+-- POWER sense -- V and I from SBAT (BATT1), carrier module supply (BATT2)
+-- IO PWM 1-2 -- flaperons via wing connectors
+-- AUX PWM -- SPUSH throttle, choke. Relay outputs -> SBAT SSR enable, SPUSH kill
+-- GPIO -- SPUSH tach in, NeoPixel data
```

Board-to-board interfaces, which is the list the connector decision is really about:

| From | To | Carries | Class |
|---|---|---|---|
| SBAT | SMAIN | HV feed 16 AWG, V sense, I sense, I2C, SSR enable, precharge end, GND | Power + signal, one harness |
| SBAT | ESC x4 | HV 12 AWG pair each | Power only, ESCs adjacent |
| SBAT | STAIL | HV spur 18 AWG | Power, runs the tail boom with the CAN trunk |
| SBAT | SPUSH | HV pair, lugs or busbar | Power |
| SMAIN | ESC x4 | CAN H, CAN L, PWM, GND | Signal, Caribou CMAIN pinout |
| SMAIN and SBAT | STAIL | One Y harness into one PTSM 4 way connector on STAIL: HV pair from SBAT, CAN pair from SMAIN, two twisted pairs in one jacket, 20 AWG | Power plus signal, decided September 28 |
| SMAIN | wing x2 | PWM, 8.4 V, GND | Mini-Fit Jr, detachable panel |
| SMAIN | SPUSH | PWM x2, relay out, tach in, 5 V, GND | Signal |
| SMAIN | shelf peripherals | GH cables to F9P, radios, altimeter, pitot | Signal, stays on the shelf |

Comparison with PT1 as flown: PT1 had one HV board (the FCHUB) and everything else as loose modules with no sensing, no tail, no servos. PT2 puts the same HV path on SBAT with the kill and sensing restored, and collapses the six loose shelf modules into SMAIN.

### 4.2 Battery decision (closed September 22)

**Pack: a scaled Longshot, 14S3P.** Decided on the September 22 battery call (Alperen, Julius, Erick, KBM, Zeynep, Thomas). Sizing was done line by line from the full 14S9P Longshot divided by about 3.2.

| Item | Value | Source |
|---|---|---|
| Configuration | 14S3P, 6.5 Ah lithium ion cylindrical cells | Sept 22 call |
| Capacity | 19.5 Ah at 14S, about 22.8 Ah 12S equivalent, roughly today's energy | Sept 22 call |
| Voltage | 58.8 V full, 42 V at 3.0 V per cell, 14S nominal 50.4 V | Li-ion cell limits |
| Mass | About 3.9 kg estimated, called 4 kg. Desk estimate with buffers | Sept 22 call |
| Connector | AS150U, which carries 4 signal pins. QS8 rejected as too hard to handle in the field | Sept 22 call |
| BMS | None at first. Iteration 1 is a dummy battery PCB with a balance lead output for the Quiver Tattu charger. Iteration 2 fits a BMS once Julius's contractor delivers one | Sept 22 call |
| Readiness | Julius: yes for the November US build if work starts the week of September 28 | Sept 22 call |
| Cost | About the same as today's pack, cells about $250 | Sept 22 call |

Why not the alternatives: Tattu energy density is about half the PT1 semi-solid pack, so equal energy in a Tattu costs about 3.3 kg, which is the whole gasoline budget (September 21 call). The Caribou 18S pack exceeds the 64 V ESC limit. Both closed.

**14S is confirmed.** Erick confirmed on the call that the MAD and KDE ESCs and the 150 KV motors run at 14S and that everything else is regulated. The 12S versus 14S research item from September 21 is therefore closed. Every voltage rating in this note now reads 58.8 V full charge instead of 50.4 V. The MAD AMPX 64 V limit leaves 5.2 V of margin.

**What SBAT carries because of this decision** (Erick's open items from September 22):

| Provision | Detail | Status |
|---|---|---|
| AS150U battery mate | Power pair plus 4 signal pins. Assign the 4 pins (CAN H, CAN L, pack temperature, pack ID or spare) with Julius so the pack PCB and SBAT agree | Required |
| BMS daughterboard header | Small header mounted board, swappable, like the X-Lite in the HESA battery boxes and the Feather approach. Connections expected: power, CAN, maybe temperature. Populated in iteration 2 | Required as a provision |
| Charging connector | Footprint for the future starter generator charge input. Not functional in November. Charging a dummy pack from a generator without a BMS is not recommended, and there is no engine in November anyway | Footprint only |
| Balance lead path | Iteration 1 charges through the main connector with a balance lead output on the pack, plugged into the Quiver Tattu charger. Nothing on SBAT unless the balance leads route through it | Confirm with Julius |
| Generator cut off | Disconnecting the generator input when the pack is full goes into the shared Longshot BMS spec, not into SBAT | Out of scope |
| Vehicle current and voltage sense | Stays on SBAT for iteration 1, because the dummy pack reports nothing. Once the BMS reports over CAN, the analog sense becomes a backup | Required |

The switching stage (fuse, precharge, MOSFET kill, SSR gating, LEDs, temperature) is unchanged by the pack choice. Re-rate it for 58.8 V. The Quiver Front parts already ran at 60.9 V on Quiver, so the precedent holds.

**Also from September 21:** Alperen's onboard BMS and in-place charging idea (balance leads to a simple BMS on the battery PCB, AS150U charge socket) is recorded as would love to have and folds into the BMS daughterboard provision above. It is dropped for this version if it threatens the November timeline.

How Quiver and Caribou mount the flight controller (verified from the board files, September 21): the FC PCB is a mezzanine daughterboard, not a side-by-side board. The Pix32 V6 FMU module plugs into J1 (AXK5S00347YG) on top of the FC PCB. Two AXK5S00047YG plugs on the bottom of the FC PCB, 52.5 mm apart, plug down into two AXK6S00547YG sockets (J3, J4) on the main board. Quiver Main V1.2 and Caribou CMAIN carry J3 and J4 at identical coordinates, so the same FC PCB fits both. CMAIN adds four SMD M3 standoffs (SMTSO3050CTJ) in a 52 x 37.5 mm rectangle around the socket pair to hold the stack down.

Open question for the call: reuse that stack unchanged on SMAIN (copy the J3 and J4 socket block and the four standoffs from CMAIN, build the existing Quiver FC PCB as-is) or merge the FC PCB circuitry into SMAIN. Recommendation is now **reuse the stack**. It costs zero FC redesign, spare Quiver FC PCBs work on Spearhead, and the 200 x 134 mm daughterboard sits above SMAIN instead of beside it. The cost is stack height against the ~90 mm shelf, and an FC PCB outline that may overhang a smaller SMAIN. Merge only if the shelf height or outline forces it.

---

## 5. Feature reuse matrix

Legend: **Copy** = take the schematic sheet and routed block as-is. **Adapt** = same circuit, changed ratings or connectors. **New** = no source.

### 5.1 SBAT, battery board

| Feature | Source | Action | Adaptation |
|---|---|---|---|
| Battery mate | Quiver Front J9, J10 custom Molex 46437 for Tattu | Adapt | Longshot 14S3P pack on AS150U (section 4.2). Replace the Molex mate with the AS150U panel footprint or 8 AWG lug landings plus a 4 pin signal header. Battery CAN to the BMS daughterboard |
| BMS daughterboard header | Feather and HESA X-Lite pattern | New | Header mounted, swappable. Power, CAN, temperature. Unpopulated in November |
| Charge connector | None | New, footprint only | Future starter generator input. Not wired in November |
| Bolt-in main fuse | Quiver Front F4 (3425L200 family), Caribou CBC AMXL-250 | Adapt | Spearhead hover 98 to 124 A, ESC-capped peak 320 A. Size the fuse against the 320 A peak, not hover. Verify the Quiver fuse rating before reuse |
| Main MOSFET switch | Quiver Front Q3 to Q7, 5x S400N08TL in parallel | Adapt | Verify the paralleled Rdson and thermal against 124 A continuous and 320 A for seconds. Part number suggests 80 V class, adequate for 50.4 V max. Add copper or count if marginal |
| Precharge and bypass | Quiver Front Q1, K-SSR, R15 and R16 (33 ohm after V1.2), Q8 bypass | Copy | Bypass moved onto the battery board in V1.2, so this block is self-contained now |
| SSR gating from FC | Quiver Front U3 CPC1019N, SSR_S from FC IO channel | Copy | Maps to the Pix32 IO channel chosen on SMAIN |
| Precharge pushbutton with LED | Quiver Front J1 | Copy | Reachable through the battery bay hatch |
| Board temperature | Quiver Front IC1 TSYS01 over I2C | Copy | Also gives the bay temperature the team wanted in July |
| Vehicle current sense | Not on Quiver Front (Quiver reads the smart battery over CAN) | New | The iteration 1 Longshot pack has no BMS. Add a shunt or hall sensor on SBAT feeding an ADC on the Pix32 power sense pins. SPH-E-001 section 4.13 gives the ArduPilot parameters |
| Pack voltage sense | FCHUB 1K:20K divider (SPH-E-001 section 4.13) | New | Resistor divider on SBAT to the Pix32 voltage sense pin |
| HV outputs | Quiver Front J11 AMASS terminal, Quiver Main J44 and J45 M6 power taps | Adapt | 4x ESC outputs at 80 A each plus one pusher HV output. XT90 or AMT0450003 screw terminals. XT60W 3-pin for ESC to motor stays in the harness |
| Status LEDs | Quiver Front LED1 to LED6 | Copy | |
| TVS on inputs | Quiver Front D1 to D9 (SMBJ13CA, SMBJ64CA) | Adapt | SMBJ64CA stands off 64 V against 58.8 V full charge, 5.2 V margin. Check the clamp voltage against the 80 V MOSFET rating, or move up one step |

### 5.2 SMAIN, cockpit board

| Feature | Source | Action | Adaptation |
|---|---|---|---|
| Pix32 V6 mounting | Quiver FC PCB stacked on the main board. CMAIN J3, J4 AXK6S00547YG sockets plus four SMTSO3050CTJ standoffs | Copy | Copy the socket and standoff block from CMAIN at the same relative geometry. Build the Quiver FC PCB unchanged. This is the single largest time saver on the project |
| 5 V and 12 V regulation | Quiver Main PS1 REC20K 5 V and PS2 REC30K 12 V, fused, with LEDs | Adapt | Caribou CMAIN uses SY8205 bucks instead. Pick one family. Caribou parts are newer and non-isolated, fine for Spearhead |
| 8.4 V servo rail, 15 A | None. PT1 used a robocombo buck | New | Adjustable synchronous buck at 8.4 V. The Kingmax CLS3015S rail for both wing servos. Design from the SY8205 block or a higher-current controller |
| CAN1 fan-out | Quiver Main J5, J7, J9, J15 and others, SM04B-GHS CAN ports | Copy | Four ESC stubs plus one tail trunk. Add a termination jumper per port end |
| CAN2 with termination switch | Quiver Main S2 and R14 (added in V1.2) | Copy | Reserved bus per SPH-E-001 section 4.9 |
| Pi 5 header, SPI, CAN for the Pi | Quiver Main J1, U1 TJA1049, U2 MCP2515, X1 | Copy | Decided September 15: Pi stays on the board |
| Ethernet switch modules and PPP-to-Ethernet ports | Quiver Main J47, J51 Gigablox, J35, J36, J41 | Copy | Needed for Starlink even without a Pi (September 15). One Gigablox may be enough, the second was added in V1.2 for payload ports |
| GPS ports | Quiver Main J7, J8 F9P DroneCAN ports, J9 and J11 M9N | Adapt | Primary F9P is serial on GPS1 for Spearhead, secondary Here4 is on CAN in the tail. Keep one serial GPS port and drop the M9N footprint |
| Telemetry and RC | Quiver Main TEL1, TEL2, TEL3 GH ports, SBUS input | Copy | SiK V3 on TELEM1, R9DS SBUS on RC IN, US-D1 altimeter on a UART, mLRS header optional |
| USB-C to FC | Quiver Main J33 | Copy | |
| Buzzer | Quiver Main LS1 | Copy | |
| Board temperature | Quiver Main IC1 TSYS01 | Copy | Second bay sensor |
| Switched 12 V payload output | Quiver Main U5 SSR, F7, F8, 12VSW | Copy | Phase 4 payload bay per SPH-E-001 section 7.1 |
| 12 V relay for payload actuator | Quiver Main K1 CPC1019N | Copy or drop | Spearhead has no dispenser. Keep the footprint as a spare FC-controlled relay for the gasoline kill wire |
| Obstacle avoidance ports | Quiver Main J5 360 LiDAR, J49 forward radar (V1.2) | Copy | Phase 3 per SPH-E-001 section 6.2 |
| Waterproof connectors | Caribou CMAIN JST JWPF family, AT13 for high current | Adapt | Spearhead is not maritime yet. JWPF is a good default for anything that leaves the shelf |
| Redundant 12 V combining | Caribou CMAIN LTC4370 x3 | Do not copy | Spearhead has one pack, nothing to combine |
| Net classes and routing rules | Caribou CMAIN `docs/routing-rules.md` | Copy | Adopt the same net-class table so DRC is consistent across Arrow boards |
| Wing connectors | None | New | 2x Molex Mini-Fit Jr 6-pin per SPH-E-001 section 4.6: PWM, 8.4 V, GND per side |
| LV kill input | None on the boards (Quiver kills at the battery board) | New | Switch input that drops the FC 5 V and peripheral 5 V rails, per SPH-E-001 section 4.3 |
| LED strip feed and GPIO | None | New | Separate 5 V feed with its own fuse, one GPIO for NeoPixel data (SPH-E-001 section 4.14) |
| KDE PWM fallback | None | New | Four ESC PWM signal pins on the ESC connectors so the KDE boards can fly PT2 if MAD never clears customs. Caribou CMAIN already puts CAN H, CAN L and PWM on one ESC connector. Copy that pinout |

### 5.3 STAIL, tail board

**Schematic complete September 28, 2026.** 81 parts, ERC 0, every part with MPN and LCSC number, BOM exported. Details in `engineering/electronics/pcbs/STAIL_PCB/README.md` and its `docs/design-log.md`. Layout waits on the tail bay envelope from Alperen and the Matek standoff height.

| Block | Decision | Source |
|---|---|---|
| Harness input J1 | One Phoenix PTSM 0,5/4 push in connector: HV+, GND, CAN_H, CAN_L. One Y harness from the fuselage, HV pair from SBAT, CAN pair from SMAIN, two twisted pairs in one jacket, 20 AWG | Quiver Main PTSM family |
| Input conditioning | 5 A nano fuse, 2x 47 uF 100 V electrolytic (ESR damps the 1.8 m cable ring), SMBJ64A unidirectional TVS (reverse polarity opens the fuse) | Quiver Front and CMAIN parts |
| 5 V logic | MP9486GN-Z 100 V buck, 1 A, fused | Caribou CBC stage, copied |
| 8 V servo rails | Two independent rails, one per ruddervator: URB4812YMD-40WR3 isolated module to 12 V, then SY8205FCC at 8.04 V (124k over 10k). PM12S-3 dropped (53 g, about $65, 60 V input limit). MP4583 single stage design kept as a backup file | Mornsun family plus CMAIN buck block |
| Matek CAN-L4-PWM | Mounted as a hat on M3 standoffs. One footprint with 34 sockets and both holes, geometry from the vendor STEP file. Bus end termination on the Matek 120R jumper, Here4 has none documented | Vendor STEP |
| Here4 | JST GH on a stub, off board | Quiver pattern |
| Servos | JR 3 pin headers, signal, 8 V, GND | |
| Status and test | Three green LEDs, twelve test points, four M3 holes | |

Rejected on the way: Tattu smart pack (mass), Mean Well modules (none in JLC stock, 53 V limit), MP9486A (1 A continuous), MAX40203 OR (single 5 V source), Micro-Fit (crimp tool).

### 5.4 Pusher interface (was SPUSH, pusher hat)

**Changed September 25, 2026 (call, Alperen, Erick, Zeynep).** There is no starter generator on PT2. The first PT2 flies an electric pusher, a gasoline version may follow as PT2.5. For the electric configuration the gas tank bay holds a secondary Longshot battery wired directly to the pusher ESC, so the pusher HV never passes through SBAT or SMAIN, the main pack does not grow, and no second Longshot size is needed. The swappable hat is dropped because the PT2 skin is multiple panels sealed with silicone, so every port has to exist from the start. SMAIN carries a permanent port set for both variants:

| Port | Variant | Notes |
|---|---|---|
| Pusher ESC port: CAN1 pair, PWM, GND | Electric | Same pattern as the four lift ESC ports. The pusher ESC reports the secondary battery voltage and current over CAN, so Alperen's request for a separate battery telemetry CAN port was not added (the Longshot has no BMS) |
| Throttle PWM | Gasoline | Quiver J24 already carries two FMU channels and 5 V |
| Engine temperature sensor, I2C | Gasoline | Alperen believes the sensor is I2C. Erick to verify against SPH-E-001 section 5.7, which assumed a MAX31855 thermocouple amplifier over SPI |
| Ignition power output, low current, XT30 class | Gasoline | 4.8 to 8.4 V per the DLE manual, switched by an FC relay output |
| Kill relay drive and tach input | Gasoline | From the September 15 call, kept as the "other requirements" Erick was to check |

The rest of this section is the gasoline interface research from September 18 and still applies to the port set above.

Gasoline variant. Engine is the proposed **DLE-55RA** (dlengine.com page and the Hobbico DLEG0455 operator manual, read September 18). Interface facts from those two sources:

| Item | DLE-55RA fact | Source |
|---|---|---|
| Engine | 55.6 cc, 5.5 HP at 8,500 rpm, idle 1,500 rpm, 7.6:1, rear carburetor with manual choke | Manual p.1 |
| Mass | Engine 1,330 g, muffler 143 g, ignition module 153 g, standoffs 128 g. About 1.75 kg before firewall and fuel | Manual p.1 |
| Thrust | 14.5 kg at 100 m altitude, 12.8 kg at 1,800 m (Ankara sits near 900 m) | dlengine.com |
| Props | 22x8, 22x10, 23x8, 23x10 | Manual p.1 |
| Ignition | Electronic CDI module. Battery lead, hall pick-up sensor lead from the engine, ignition wire to an NGK CM6 plug at 0.45 to 0.51 mm gap | Manual p.1, p.7 |
| Ignition supply | "Any 4.8-8.4V, 1000mAh and above capacity battery." SPH-E-001 section 5.4 says 5 to 6 V, correct it to 4.8 to 8.4 V | Manual p.6 |
| Ignition switch | Red 3-pin connector lead for a manual on/off switch. The manual says always install one, accessible from outside the aircraft | Manual p.2, p.4, p.6 |
| Remote kill | Optional opto-isolated TX kill switch (DLEG9205) wired between the manual switch and the ignition, with an external red LED that shows the ignition is armed. Both the manual switch and the kill must be on before ignition can occur | Manual p.4, p.6, p.7 |
| RPM output | The ignition module ships "w/ additional tachometer lead", labelled "Tachometer Lead / RPM Signal Output" on the wiring diagram. Signal level and pulses per rev are not stated. Measure on the bench before assigning a Pix32 pin | Manual p.2, p.7 |
| RF and placement | Install all electronic radio devices at least 305 mm from the engine. Throttle and choke servos at least 305 mm away, non-metallic pushrods. Do not put the throttle servo in the engine compartment | Manual p.2, p.7, p.9 |
| Firewall | At least 9.5 mm 5-ply plywood, glued and doweled. Engine length from standoff face to prop washer is 171 mm | Manual p.7 |
| Starter | "A 24 volt electric starter is recommended" for hand starting with a helper holding the aircraft. This is a handheld field starter, not an onboard unit. No onboard starter exists for this engine | Manual p.10 |
| Fuel | 87 to 93 octane at 30:1. SPH-E-001 section 5.1 says 40:1 typical, correct it | Manual p.1 |

What SPUSH therefore carries for the gasoline variant:

| Feature | Source | Action | Adaptation |
|---|---|---|---|
| Ignition supply 4.8 to 8.4 V, isolated from avionics | Quiver Main REC isolated converter block | Adapt | Isolation matters because the CDI pulses on its supply. Fuse it. The manual's own answer is a separate 1,000 mAh pack, which is also acceptable for PT2 |
| Manual ignition switch input | None | New | Panel switch on the fuselage, in series before everything else, per the manual |
| Remote kill | Quiver Main K1 CPC1019N relay block, or the DLEG9205 opto kill as-is | Copy or buy | Wired between the manual switch and the CDI. Driven by a Pix32 relay output. Add the armed LED to the exterior |
| Tach input | None | New | One protected GPIO with pull-up, RPM1 in ArduPilot. Bench-characterise the module's tach lead first |
| Throttle and choke servo PWM | Tail servo header pattern | Copy | Two JR headers, plain PWM from Pix32 aux outputs, servos mounted 305 mm or more from the engine |
| Starter | Nothing onboard | Drop | Remove the starter relay and 2S battery footprints from the earlier plan. Field starter only |
| CHT and EGT | None | New footprint only | MAX31855 over SPI, SPH-E-001 section 5.7 |
| Electric variant | HV XT90 out plus ESC signal and telemetry | New | Same outline, different population |

rusEFI (open source ECU, raised as context): it is a fuel-injection engine management platform for automotive-scale engines. The DLE-55RA is carbureted with its own CDI, so rusEFI has nothing to control here. Not applicable unless a future engine is fuel injected.

### 5.5 Future items from Thomas's September 6 notes

Each item mapped to a board so nothing gets lost. Phase per SPH-E-001 unless noted.

| Item | Board | Phase | Source block | Note |
|---|---|---|---|---|
| Second GPS plus CAN node driving tail servos | STAIL | PT2 | Section 5.3 | Already scoped |
| Aileron (flaperon) servos wired into the cockpit | SMAIN | PT2 | Wing connectors, 8.4 V rail | Already scoped |
| Electric pusher and ESC | SPUSH electric variant | PT2 | Section 5.4 | Already scoped |
| Gas engine pusher | SPUSH gasoline variant | PT2 to PT3 | Section 5.4 | The interface table above answers "what does this need" |
| Arrow-style Pix32 mount | SMAIN | PT2 | Quiver FC carrier | Already scoped |
| Precharge and battery current monitor | SBAT | PT2 | Quiver Front PCB | Already scoped |
| SIYI air unit (HM30) and SIYI camera (A8) | SMAIN | PT3 | Quiver Main HM30 telemetry ports and J2 Ethernet to the A8 gimbal | Copy the Quiver ports. Adds the video link Zeynep asked for after STORK |
| Payload interface, Quiver style plus HV supply | SMAIN and SBAT | Phase 4 | Quiver PT3 spec 10-pin payload connector (4 Ethernet, 2 CAN, 1 analog, 12 V and GND, 1 spare) plus a fused HV output on SBAT | Aligns with SPH-E-001 section 7 |
| Ethernet | SMAIN | PT2 | Quiver Main Gigablox | Already scoped |
| Nav and position lights, beacon, strobes | SMAIN | PT2 | LED feed and NeoPixel GPIO (SPH-E-001 section 4.14) | Thomas asked "CAN?". A DroneCAN LED node would put the strips at the wingtips with only CAN and power on the wing connector. Decide when the wing connector pinout is fixed |
| Remote ID | SMAIN | PT3 | None | ArduPilot OpenDroneID module on CAN or serial. Reserve one CAN port and a 5 V feed. Regulatory need depends on where PT3 flies |

---

## 6. Size and placement constraints

| Item | Number | Source | Status |
|---|---|---|---|
| Avionics shelf | ~180 x 250 mm footprint, ~90 mm tall | June 12 estimate, SPH-E-001 section 2.3 | PT2 fuselage is being redesigned. Get new numbers from Alperen before fixing any outline |
| PT1 fuselage exterior | 240 mm wide by 275 mm tall, constant over a 540 mm mid body. Nose 355 mm, aft 200 mm | Aerodynamics information note 0001 and `fuselage.step`, May 25, 2026 | Outer mold line, so usable interior width is under 240 mm. PT2 fuselage is a new design, so this is a sanity bound and not a requirement |
| Quiver Main | 356 x 254 mm | File | About twice the shelf area. SMAIN must shrink to roughly half. Its 254 mm short side is wider than the whole PT1 fuselage |
| Quiver FC carrier | 200 x 134 mm | File | Fits the shelf alone, not alongside a Quiver-size main |
| Quiver Front | 248 x 168 mm | File | Larger than the PT1 pack face (190 x 78 mm). The Longshot 14S3P envelope is not yet drawn, Julius starts it the week of September 28. SBAT shrinks once the Molex mate goes, and grows back a little for the BMS header and charge connector |
| Tail bay | Unknown | Alperen | NACA 0018 today, may go to 0021 for depth |
| ESC location | In the fuselage, accepted September 15 by Alperen | Call | Confirm with the PT1 build photos. Sets whether SBAT carries the ESC outputs and where the heat goes |

---

## 7. Decisions to close on the September 18 call

1. **Board count and FC mounting.** SBAT, SMAIN, STAIL, SPUSH, plus the existing Quiver FC PCB stacked on SMAIN through the J3 and J4 mezzanine sockets, as Quiver and Caribou do. Recommendation: reuse the stack unchanged. Merge the FC circuitry into SMAIN only if shelf height or outline forces it.
2. **Hardware HV kill returns on PT2** via the SBAT MOSFET stage. Confirm this closes the section 9.3 deviation and that Zeynep is fine with the precharge and SSR sequence on the Pix32.
3. **ESC location** confirmed as fuselage. Then SBAT carries four 80 A outputs and the ESCs mount next to it.
4. **Regulator family.** Quiver REC isolated modules or Caribou SY8205 bucks. One family across all Spearhead boards.
5. **Connector family for off-shelf runs.** JST JWPF (Caribou) versus JST GH (Quiver). Recommendation: JWPF for anything that leaves the shelf, GH on the shelf.
6. **Pi and Ethernet stay** (already agreed September 15). Confirm one or two Gigablox modules.
7. **Starlink weight.** SPH-E-001 section 6.1 rejected Starlink Mini at about 1 kg. Thomas's Pi-for-Starlink plan needs that reconciled or the Pi is carrying a dead use case.
8. **Starter-generator stays out** (agreed September 15). Footprint-only provisions on SPUSH are acceptable.
9. **Battery.** Closed September 22: scaled Longshot 14S3P, about 4 kg, AS150U, BMS daughterboard provision on SBAT (section 4.2). Design every board for 58.8 V full charge.
10. **Dimensions.** Alperen delivers shelf, battery bay and tail bay envelopes for the PT2 fuselage.
11. **Lead time.** Julius confirms fab plus assembly plus shipping to the US.

---

## 8. Inputs still needed

| From | Item |
|---|---|
| Alperen | PT2 shelf, battery bay and tail bay dimensions. PT1 build photos of the bay. Which ESC set flies PT2 (MAD or KDE) |
| Zeynep | Pix32 pin plan for the gasoline engine (tach, kill, throttle, choke), ArduPilot ICE parameters, the CAN-L4-PWM servo rail continuity test |
| Thomas | Notebook notes from the Ankara visit on what the PCB needs |
| Julius | Fab and assembly lead time, and confirmation of which Quiver Main stackup was fabricated (4 or 6 layers) |
| Erick | Current sense method on SBAT (shunt vs hall), 8.4 V rail design, fuse and MOSFET verification against 320 A, bench characterisation of the DLE tach lead |
| Alperen | Which DLE variant is on order (55RA per Thomas), and whether the DLEG9205 opto kill and a tachometer come with it |

---

## 9. Working method for the reuse

- Open every source board in KiCad 10. Save copies into the Spearhead repo under `engineering/electronics/pcbs/<BOARD>/kicad/` using the Caribou folder convention, with symbols, footprints and 3D models copied into `kicad/libs/` so nothing depends on a personal library.
- Schematic reuse: the Quiver Main project already has hierarchical sheets (`CAN circuit.kicad_sch`, `ETHERNET.kicad_sch`). Import those sheets whole. Add the FC carrier schematic as a sheet.
- Layout reuse: select the routed block in the source board (footprints plus tracks plus vias), copy, paste into the new board. KiCad keeps the routing. Re-annotate afterwards and re-run the netlist sync.
- Adopt the Caribou net-class table before routing anything new so DRC is consistent.
- Run DRC and ERC from the command line after every session:

```
kicad-cli pcb drc --format json --output drc.json <board>.kicad_pcb
kicad-cli sch erc --format json --output erc.json <board>.kicad_sch
```

- Keep a per-board README and docs folder as Caribou does, and log every departure from the source board in a `docs/reuse-log.md` so a reviewer can trace what changed and why.

### 9.1 Design start order

Base project for SMAIN, decided September 21: **Quiver Main V1.2**, not Caribou CMAIN. CMAIN is the cleaner file set (KiCad 10, project-local libraries, ERC 0, DRC 2 cosmetic), but it takes 12 V from six Caribou battery boards and never sees battery voltage. Its power supply sheet is 108 parts of LTC4370 combiners and 30 V bucks that Spearhead cannot use from an HV feed. Quiver Main has the Spearhead power architecture: one HV feed, isolated REC30K-4812SZ, REC20K-4805SZ and URB4805S-6WR3 converters, four ESC ports, 164 parts against 248. Its file hygiene gap (KiCad 9 format, links to personal global libraries) is a one-time fix. CMAIN stays the donor for JWPF connectors, SY8205 point-of-load bucks, spare SSR channels and the net class table. The Quiver FC PCB is built unchanged and needs no design work.

| Order | Board | Method | Why this position |
|---|---|---|---|
| 1 | STAIL | New small project. One HV input stage, one 8 V buck, one 5 V rail, CAN in and out with a termination jumper, CAN-L4-PWM header, Here4 connector, two JR headers. The HV to 8 V stage is a new design, because the CMAIN SY8205FCC buck is a 30 V input, 5 A part. SY8205 is reused only for the 5 V rail, fed from the 8 V rail | About 25 parts. Exercises the full flow (schematic, footprint assignment, layout, DRC, fab outputs) on a board where a mistake is cheap. Does not depend on any open decision |
| 2 | SMAIN | Copy of the Quiver Main V1.2 project, renamed and upgraded to KiCad 10 format. Localise the libraries once, then edit sheet by sheet. Lock J3 and J4 before Update PCB from Schematic. Decide six or four copper layers before rerouting. Shrink the outline last | Same power architecture as Spearhead, so the power sheets stay. New circuits: 8.4 V servo rail, wing connectors, SPUSH signal connector, SBAT harness connector. Imports from CMAIN: the four FC hold-down standoffs around J3 and J4 (Quiver Main has none, its four standoffs are the Pi 5 mount), JWPF connectors, SY8205 point-of-load bucks, net class table |
| 3 | SPUSH | New small project, two population variants on one outline | Small. Gasoline side waits on the DLE tach lead bench measurement |
| 4 | SBAT | Copy the Quiver Front project, migrate it to KiCad 10 and project-local libraries, swap the battery mate, add the sense block | Depends on the September 22 battery call. Switching stage first, battery mate and sense block last |

Scaffold, September 21: the four projects exist under `engineering/electronics/pcbs/` with a README per board, reuse logs for the two copied boards and baseline ERC and DRC reports. SMAIN matches Quiver Main V1.2 and SBAT matches Quiver Front by net list comparison (164 and 87 components, identical reference designators and nets). Both were upgraded to KiCad 10 format and still link to the original designer's global libraries, which each README explains how to localise. The inherited net classes stop at Power_5A with 0.3 mm clearance, so SBAT needs an HV class (IPC-2221, 0.6 mm from 31 to 100 V, uncoated external).

Per-board gate before layout starts: schematic passes ERC clean, every symbol has a footprint from `libs/`, connector pinouts are written in the board README, and a second person has read the schematic PDF.

---

## 10. Revision history

| Rev | Date | Author | Changes |
|---|---|---|---|
| 0.1 | September 18, 2026 | errrks | Initial draft for the September 18 call. Board inventory read from the KiCad files on this date. Scope from the September 15 call. Added Thomas's September 6 as-built and future lists, and the DLE-55RA interface facts from the dlengine.com page and the DLEG0455 manual |
| 0.2 | September 21, 2026 | errrks | FC mounting verified as a stacked mezzanine, recommendation changed to reuse the stack. Tattu smart pack considered and dropped, Longshot pack questions added (section 4.2). Design start order added (section 9.1) with the CMAIN ERC and DRC baseline. KiCad projects scaffolded. Corrected two reuse claims: SY8205 cannot take the HV feed, and CMAIN has no HV front end, which moved the SMAIN base to Quiver Main V1.2 |
| 0.3 | September 25, 2026 | errrks | Folded in the September 21 and 22 calls from the wiki. Battery closed as a scaled Longshot 14S3P on AS150U with a BMS daughterboard provision, 12S replaced by 14S throughout, SBAT rows for the BMS header and charge connector, TVS margin at 58.8 V. Board set confirmed on the September 21 call (FC, main, battery with motor and tail outputs, optional pusher, tail) |
| 0.4 | September 28, 2026 | errrks | STAIL schematic complete, section 5.3 rewritten as decisions. Tail harness is one PTSM 4 way connector fed by a Y harness (section 4.1 interface table updated). Schematic conventions and the check script recorded in `engineering/electronics/pcbs/README.md` and `engineering/electronics/tools/sch_check.py` |
| 0.5 | September 28, 2026 | errrks | Folded in the September 25 call: SPUSH hat dropped, electric pusher on a secondary Longshot battery wired directly to the pusher ESC, permanent pusher port set on SMAIN (section 5.4, board set table). SMAIN pin plan Rev 1 written and decided (`SMAIN_PCB/docs/pinout.md`). Section 4.1 diagrams and the interface table still show the SPUSH hat and the SBAT pusher output, update on the next pass |
