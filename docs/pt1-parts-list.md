# Spearhead PT1 Electrical Parts List

**Doc ID:** SPH-E-003
**Author:** errrks
**Date:** June 17, 2026
**Status:** Ready for order review (next Spearhead call)
**Reference:** SPH-E-001 (electrical master, Rev 0.18) §4.16 cost baseline, §4.6 harness connection table

Procurement list for the PT1 order cycle, per the June 11 call. Focus areas requested by Alperen: voltage regulators, the tail CAN-to-PWM board, wiring and connectors, and items with long lead times.

**Sourcing rule:** buy in Turkey where stock exists (battery at motorobit, FC at rx-dynamic, servos at thkmodelucak, tail CAN board at Aykut Havacılık), then EU stock (Matek EU resellers, Holybro EU, Aeroboticshop NL) to keep customs simple, same as the MAD order with the ATR certificate. The Ainstein US-D1 is believed to be on hand in Turkey already, confirm before any US order.

---

## Priority 1 — order or request quote this week (long lead or gating)

| # | Item | Spec / model | Qty | Source candidate | Est. unit | Lead risk | Notes |
|---|---|---|---|---|---|---|---|
| 1 | [Flight controller](https://rx-dynamic.com/urun/holybro-pixhawk-6c-pm02/) | **Pixhawk 6C + PM02 V3** (full-size 6C, not Mini; analog PM02 V3 12S power module + cable set, **no GPS in this bundle**) | 1 | rx-dynamic (TR), ₺19,900 (~$600) | ~$600 | Medium | **D4 reopened and reclosed June 17: full 6C selected** (supersedes the June 11 6X). Its IO PWM still drives the flaperons (wing CAN node stays dropped), it has 5 UARTs (vs the Mini's 4) for the serial F9P, and its **two analog POWER ports** let the FCHUB vehicle V/I feed POWER2 directly (SPH-E-001 §4.13). No M9N in this bundle (a +M9NGPS 6C variant exists separately) — the primary GPS is the F9P (line 4b) |
| 2 | [CAN-to-PWM](https://www.aykuthavacilik.com/urun/ap-periph-can-node-can-l4-pwm) | **Matek CAN-L4-PWM** (AP_Periph, 9 PWM, JST-GH CAN, 3.5g) | 2 | Aykut Havacılık (TR), 2,060 TL (~$50) | $50 | Low (TR stock) | 1 install (tail, node 12, ruddervators) + 1 spare. The wing unit was dropped, the FC IO PWM drives the flaperons. Gates WP-E04 bench validation. Alternate at the same shop: CAN-L431 (5 PWM, 3 UART, I2C, dual parallel GH-4P CAN, 2,500 TL) if the L4-PWM sells out |
| 3 | [Tail servo regulator](https://www.aykuthavacilik.com/urun/power-module-12s-w-3xbec) | **Matek PM12S-3** (9–55V in w/ TVS. Vx selectable 5.25/6/8V at 15A cont / 25A peak, fixed 5.2V/4A, fixed 12V/4A, 1K:20K divider, 53g) | 1 | Aykut Havacılık (TR), 6,900 TL (~$170) | $170 | Low (TR stock) | **Tail only now (June 17):** Vx at 8V for the ruddervator Kingmax, fixed 5V for the tail adapter + Here4. Kept because it supplies both the 8V servo rail and the 5V logic rail in one TVS module. The forward unit is dropped — the wing rail moves to a robocombo buck (line 3b). **Price flag:** TR price ~4× Matek direct, accepted for customs simplicity |
| 3b | [Wing servo regulator](https://www.robocombo.com/dc-dc-8-60v-15a-ayarlanabilir-voltaj-dusurucu-modu-4017) | **robocombo DC-DC buck** (8–60V in, 15A adjustable) | 2 | robocombo (TR) | ~$15 | Low | 1 install (wing flaperon rail, HV input at the FCHUB pads, output set to **8.4V** for full Kingmax torque) + 1 spare. 8.4V is why this replaces the forward PM12S-3 (Vx caps at 8V). Lock the output pot, confirm true synchronous 15A with heatsink, scope for overshoot before connecting a servo |
| 4 | Tail GPS | **CubePilot Here4** (secondary GPS, CAN node 11) | 1 | Aeroboticshop (EU) | $300 | Medium | Secondary GPS. D1 closed June 17 with a serial F9P primary, so this stays a single Here4 (no matched pair) |
| 4b | Front GPS | **u-blox F9P** on a serial UART (GPS1), primary | 1 | Salvaged F9P on hand; antennas on hand | $0 | — | **D1 closed June 17:** salvaged F9P on serial (compass yaw), not DroneCAN. Uses the on-hand antennas. No purchase needed |
| 5 | Radar altimeter | **Ainstein US-D1** (UART + CAN, 100 Hz, 50 m, 110g, 5–13V) | 1 | **Believed on hand in Turkey, confirm unit + interface variant** | $0 if found | On-hand check | Selected June 11 over TF03-180 (radar tolerates prop wash dust). Only if the unit does not turn up: quote via Ainstein/distributor, then it becomes the longest lead item (~$280, High) |
| 6 | [Battery](https://www.motorobit.com/444v-12s-22000mah-15c-solid-state-lipo-batarya) | **motorobit 12S 22,000 mAh 15C solid-state LiPo** (44.4V, 3,709g, 190×78×126mm, QS8-S socket) | 1 | motorobit.com (TR) | TBD | Low (TR) | **D5 reopened/reclosed June 17:** single solid-state pack, supersedes the ProFuse 16Ah. 3,709g is ~113g lighter than the ProFuse so no CG penalty for the +6Ah. Send 190×78×126mm to Alperen for the forward bay. 15C = 330A continuous (covers the 320A peak, more sag than 60C). **QS8-S is anti-spark** (resolves the §4.3 connector question). Confirm price |
| 7 | [PDB](https://www.aykuthavacilik.com/urun/fchub-12s-v2-w-curr440a-5v-12v) | **Matek FCHUB-12S V2** (8–60V, 4×70A cont / 4×110A burst, 5V/5A + 12V/4A + 3.3V BECs, 440A sensor, VBat 1K:20K divider, PINIO-switchable 12V/2A pad) | 2 | Matek / EU reseller | $30 | Medium (confirm EU stock) | Promoted to P1 June 11: the only HV distribution board, it gates the entire HV path and every BEC rail, so a stock surprise stalls the build. 1 install + 1 spare. V2 divider closes the voltage sense question, the 12VSW pad can drive the HV kill contactor coil (D2) |

## Priority 2 — order with the main batch

| # | Item | Spec / model | Qty | Source candidate | Est. unit | Notes |
|---|---|---|---|---|---|---|
| 7 | [Control surface servos](https://www.thkmodelucak.com/urun/kingmax-servo-cls3015s-high-torque-8-4v-1958) | **Kingmax CLS3015S** (HV 6.0–8.4V, 35 kg·cm @ 8.4V, 0.15s/60°, metal gear, dual BB, waterproof, 80g, 25T) | 6 | thkmodelucak (TR), 3,950 TL each (~$110) | ~$110 | Medium | **Moved to final HV servos on PT1 (June 17), 5V test stage dropped.** 4 install (2 flaperon @ 8.4V robocombo + 2 ruddervator @ 8V PM12S-3) + 2 spares. One SKU for all surfaces. Note +100g aft tail mass for CG |
| 8 | CAN hub | CUAV CAN HUB or equivalent passive JST-GH splitter | 1 | CUAV / EU reseller | $15 | Fans CAN1 out to the four ESC stubs + nose + tail trunk |
| 9 | CAN terminators | 120Ω JST-GH plug | 2 | with hub order | $3 | FC end + last tail node, if not built into hub/devices |
| 10 | Battery connectors | QS8-S plug ×3 (mates the pack socket) + QS8-S ×2 (charge lead and bench adapter) | 1 lot | Amass / TR hobby | $20 | The pack's QS8-S is the ground safing disconnect and is anti-spark, so no AS150U adapter is needed and there is no mating spark regardless of the D2 kill choice |
| 11 | ESC junction connectors | XT90S pairs | 8 | Amass / EU reseller | $6 | Boom exceeds the 800mm factory lead (confirmed June 11), so all four ESC power leads are extended at the boom-fuselage junction with XT90S. On the tail-carrying boom this junction also breaks out the 18 AWG tail HV feed (H5). Qty covers 4 installs + spares/bench. Exact added length pending boom dimensions from Alperen |
| 12 | LV connectors | XT60 pairs ×4, XT30 pairs ×4 | 8 | Amass / EU reseller | $3 | XT60: power module leg. XT30: BEC inputs |
| 13 | Wing panel connectors | Molex Mini-Fit Jr 6-pin sets | 3 | Mouser/Digikey | $10 | 2 install + 1 spare, pending OQ-09 confirmation (cheap, order anyway) |
| 14 | HV wire | **8 AWG silicone, 1 m red + 1 m black** (battery harness) and **12 AWG silicone, 5 m red + 5 m black** (ESC lead extensions) | 1 lot | local / EU | $45 | 8 AWG for the battery QS8-S to FCHUB pads: 4 AWG does not fit the QS8 solder cups, and on this short run the connector (QS8 ~120A class), not the wire, is the thermal limit — monitor connector temp at worst-case hover (~124A). 12 AWG matches the AMPX factory leads and is the confirmed extension gauge throughout: MAD rates the factory 12 AWG at 80A and hover draw is only ~25–31A/motor, so no upsize is needed |
| 15 | LV + signal wire | 16 AWG (power module pigtail), 18 AWG (tail HV spur, 5 m), 20–26 AWG assortment, 24 AWG twisted pair (CAN, 5 m) | 1 lot | local / EU | $40 | Check the JST-GH cable stock first, the team already has many Pixhawk cables |
| 16 | Main harness parts | 8 AWG terminations, strain relief, cable lacing | 1 lot | local | $15 | Series link dropped, single 12S pack selected |
| 17 | LV kill switch | Illuminated rocker or key switch, 10A | 1 | local | $15 | Independent of the HV kill decision (D2), order now |
| 18 | LED system | WS2812B strip + 22 AWG feed wire | 1 lot | local / EU | $20 | E-REQ-11. Own 5V feed, off the peripheral rail |
| 19 | Pitot transducer backup | Holybro / Matek I2C DLVR airspeed sensor | 1 | EU reseller | $30 | Backup in case the drawer pitot (OQ-06) has no usable transducer |
| 20 | Cell checker | Plug-in balance checker/alarm | 2 | local | $5 | Ground checks only, the dumb pack has no BMS |
| 21 | Consumables | Solder-seal sleeve assortment (LV sizes), dual wall adhesive heat-shrink assortment, regular heat-shrink assortment, braided sleeve, zip ties, mounting hardware, threadlock | 1 lot | local | $50 | Solder sleeves are for LV and signal splices. HV splices get a soldered lap joint under dual wall adhesive heat-shrink |
| 22 | 12 AWG splice sleeves | Solder-seal sleeves sized for 12 AWG (ESC lead extensions) | 12 | local / EU | $10 | On the 80A path prefer a soldered lap joint + dual wall adhesive heat-shrink, use sleeves only where iron access is awkward, and pull-test every joint |
| 23 | Motor phase bullets | 5mm gold bullet sets (mate the MAD phase leads) | 12 | local / EU | $2 | Only where phase leads get connectorized instead of soldered direct |
| 24 | JST kits | Genuine **JST-GH 1.25mm** kit with pre-crimped pigtails (CAN, Pixhawk peripherals) + **JST-PH 2.0mm** kit | 2 kits | Mouser/Digikey | $25 | Generic "1.25mm Micro JST / PicoBlade" is not GH and will not latch, buy actual GH. The FCHUB V2's JST-SH 1.0mm cable ships with the board |
| 25 | Servo leads | JR 3-pin extensions and pre-crimped pigtails, assorted lengths | 1 lot | local TR hobby | $15 | Tail bay servo runs (H15/H16), wing root runs, and bench |

### Battery interface build (how the PM02 V3 pairs to the battery)

No Y-cable. The FCHUB-12S battery pads are the split point:

1. Battery main connector (pack default) → 8 AWG red/black, ~0.4 m → solder to the FCHUB B+ / B− pads. The HV kill (D2) goes inline in this run.
2. At the same pads, solder a 16 AWG pigtail ending in an XT60 female → PM02 V3 XT60 input. The PM02 V3 sees pack voltage at the same electrical node a Y-splitter would give it.
3. PM02 V3 XT60 output: unused, insulate it. Its 6-pin CLIK-Mate cable goes to Pixhawk 6C POWER1 (5.2V analog, BATT2).
4. **Vehicle V/I:** a separate 4-wire CLIK-Mate from the FCHUB 440A sensor + VBat divider goes to 6C **POWER2** (Curr→CURRENT2, V→VOLTAGE2, GND→GND, POWER2 5V pins open). FCHUB = BATT1 (vehicle current + SoC), enabled by the 6C's second analog port. No cutting the PM02 harness. See SPH-E-001 §4.13.

Fallback if two wire sets crowd the pads: solder both the 8 AWG main and the 16 AWG branch into the main connector's solder cups together (8 AWG + 16 AWG fit a 6mm cup), which moves the split into the connector barrel. Electrically identical. Pick one approach, not both.

## Priority 3 — blocked on a decision (do not order yet)

| ID | Decision | Options on the table | Owner | Needed by |
|---|---|---|---|---|
| D1 | ~~Nose (primary) GPS, OQ-12~~ | **Closed June 17: u-blox F9P on a serial UART (GPS1)**, compass yaw, on-hand antennas (line 4b). Mixed serial-F9P + CAN-Here4 means no moving-baseline GPS yaw | Team call | Closed |
| D2 | ~~HV kill implementation, OQ-05~~ | **Closed June 17: no hardware HV kill for PT1.** Software E-stop (ArduPilot motor emergency stop on an RC switch, AMPX 2 s CAN watchdog backstop) + QS8-S anti-spark unplug for ground safing. Recorded §9.3 deviation (E-REQ-02). LV kill retained. Rationale: no clean contactor placement for PT1 without added mass, and the QS8-S anti-spark already manages connection inrush. Contactor/MOSFET analysis kept as record in SPH-E-001 §4.3; full HV kill returns for the final product | Erick + Alperen | Closed |
| D3 | RC link + telemetry radio, OQ-04 | SIYI HM30 (shared with Quiver?) vs ExpressLRS + SiK/RFD900x. Receiver and radio lines wait on this | Team call | Before harness finalization |
| D4 | FC selection | **Reclosed June 17: full-size Pixhawk 6C + PM02 V3** (line 1), supersedes the June 11 6X. IO PWM still drives the flaperons, 5 UARTs fit the serial F9P, and two analog POWER ports carry the FCHUB V/I. Sourced from rx-dynamic TR | Alperen (sourcing) | Closed |
| D5 | Battery format | **Closed June 17: single motorobit 12S 22Ah 15C solid-state pack** (line 6), 3,709g, 190×78×126mm, QS8-S anti-spark socket. Supersedes the June 11 ProFuse. Only price left to confirm | — | Closed |

## Ordered — in transit (MAD Motor Poland, no action except customs paperwork)

Ordered ~May 1, 2026 with the ATR certificate requested for zero-tax EU-to-Turkey import. Status June 11: MAD sent the purchaser/customs document, Alperen must print, sign, and return it. Listed here for completeness, full specs in SPH-E-001 §4.5.

| Item | Spec / model | Qty | Unit | Total | Notes |
|---|---|---|---|---|---|
| [VTOL motors](https://store.mad-motor.com/products/v8013-pro-ipe-feathering-propeller-autocenter-evtol-drone-motor) | **MAD V8013 PRO IPE 150KV** (591g, 15 kgf max thrust, prop autocenter indexing, 12S) | 6 | ~$185 | ~$1,110 | 4 install + 2 spares |
| [ESCs](https://store.mad-motor.com/products/ampx-80a-5-14s-drone-esc) | **MAD AMPX 80A V2 DroneCAN** (5–14S, IP67, CAN telemetry, 12AWG/800mm power leads) | 6 | $60 | $360 | 4 install + 2 spares. CAN nodes 1–4 |
| [Propellers](https://store.mad-motor.com/products/fluxer-pro-26-7-8inch-matt-carbon-fiber-propeller-for-the-long-range-flight-time-professional-hexacopter) | **MAD FLUXER PRO 26×7.8 MATT** carbon, fixed (69g/blade, Ø10mm bore, Ø20mm M3×4 hub) | 6 pairs | $193.90 | ~$1,163 + $60 ship | 12 blades (6 CW + 6 CCW), 4 installed, spares absorb early breakage. Ultralight variant deferred to ~Aug 2026 |
| **Subtotal (ordered)** | | | | **~$2,633** | |

## On hand — no order

| Item | Status |
|---|---|
| GPS antennas | On hand per June 11 call |
| Pixhawk / JST-GH cable assortment | On hand, top up only after stock check (line 15) |
| Pitot tube | In Alperen's drawer, model and transducer TBD (OQ-06) |

---

## Totals (estimate, decisions pending)

| Block | Est. |
|---|---|
| Priority 1 (6C + PM02 V3 set ~$600, single Here4, F9P TBD, motorobit pack TBD, 1× PM12S-3, 2× robocombo buck, 2× CAN-L4-PWM, 2× FCHUB-12S V2, US-D1 on hand) | ~$1,300 + battery/F9P TBD |
| Priority 2 (incl. 6× Kingmax CLS3015S ~$660) | ~$1,000 |
| Decision items (D2 contactor route, if chosen) | ~$160 |
| If the US-D1 does not turn up in Turkey | +$280 |
| **PT1 order cycle total** | **~$2,300–2,750 + battery + F9P** |
| Already ordered and in transit (MAD propulsion set, memo only) | ~$2,633 |

Cost baseline and per-phase totals stay in SPH-E-001 §4.16. Cable-level connection detail (what plugs into what, AWG, lengths) is SPH-E-001 §4.6, harness connection table.
