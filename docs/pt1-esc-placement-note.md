---
title: "Spearhead PT1 KDE ESC Placement Research Note"
sidebar_label: "PT1 ESC Placement"
sidebar_position: 8
---

# Spearhead PT1 KDE ESC Placement Research Note

Internal doc ID: SPH-E-007 (renumbered October 16, 2026, was SPH-E-005). Companion to SPH-E-001 (Spearhead Electrical Master), SPH-E-002 (PT1 architecture and decisions), and SPH-E-003 (PT1 parts list).

## Status

`Draft — pending team decision`

`Revision History: Rev 0 (July 27, 2026): initial research note answering the July 24 call action item on ESC placement for the KDE fallback propulsion set.`

`Reference: SPH-E-001, SPH-E-002, SPH-E-003, KDE Direct UAS55HVC product documentation`

## Context

Turkish customs rejected the MAD Motor shipment on July 17, 2026. The PT1 working plan (July 24 call) is to fly the KDE motors and KDE-UAS55HVC ESCs from the Spearhead testbed. This changes the ESC interface from DroneCAN to PWM and reopens the ESC placement question, because the KDE arrangement Alperen is currently designing puts each ESC beside its motor, which turns the boom run into a long DC power run. KDE places a hard 12 inch limit on DC lead length and requires added capacitance beyond it.

This note evaluates three placements and recommends one. The MAD architecture in SPH-E-001 (AMPX ESCs on booms, DroneCAN) is unchanged as the design intent and returns when the re-declared shipment clears.

## KDE constraints (vendor guidance)

From the UAS55HVC product page, the KDE electronics FAQ, and the KDE technical article "The Dangers of Long DC Wires":

| Item | Value |
|---|---|
| Voltage range | 3S to 12S LiHV (max 52.2 V) |
| Current | 55 A continuous (180 s), 95 A peak (5 s) |
| Power leads | 13 AWG silicone, bullet connectors |
| Motor leads | 14 AWG silicone, 3.5 mm / 4.0 mm bullets |
| Control lead | 22 AWG 3-wire JR, PWM 50 to 600 Hz adaptive |
| Max DC wire length | **12 inches** between ESC and power source |
| Cooling | 5 mph minimum airflow across the heatsink for full ratings, 95 °C max |
| Ingress | IP56 |

Verbatim rules that drive this decision:

- "Extending the ESC power leads is at your own risk and voids the warranty. Additional capacitance must be added if the DC wire length between the ESC and power source is longer than 12 inches."
- Capacitance sizing: **2.35 µF per amp for each inch of DC wire beyond 12 inches** (KDE's worked example sizes at the ESC rating, not at typical draw).
- Capacitors must sit within 4 inches of the ESC input caps, split across several parallel cans, low ESR, derated. If the total DC run exceeds 3 feet, capacitance must additionally be distributed every 8 to 12 inches along the wire.
- KDE's stated preference: "Keep the DC wires short and extend the AC wires between the motor and ESC." The FAQ repeats it: "We recommend extending the ESC motor wires instead."

## Spearhead geometry

From SPH-E-001 §4.14 (H3) and the harness table: the FCHUB to motor run is ~1.4 m on the front booms and ~1.6 to 1.7 m on the rear booms (rear booms run 20 to 30 cm longer). That is 55 to 67 inches, which exceeds both the 12 inch capacitance threshold and the 3 foot distributed-capacitance threshold by a wide margin. Exact lengths remain pending Alperen's final boom dimensions.

Hover current with the KDE set is unverified. The MAD sizing estimated 25 to 31 A per motor at hover. The KDE motors are smaller with a ~4.5 kgf continuous thrust ceiling (recalled, unverified), so hover throttle and current will sit higher, plausibly 30 to 40 A. A bench run on the testbed with a current measurement is the only way to pin this down.

## Option A: ESC at the motor, long DC run (current mount design)

This is the arrangement in Alperen's July 24 mount sketches.

Capacitance required, sized at the 55 A rating per KDE's method:

| Boom | Run | Length beyond 12 in | Capacitance per ESC |
|---|---|---|---|
| Front (~1.4 m / 55 in) | 55 in | 43 in | 2.35 × 55 × 43 ≈ **5,600 µF** |
| Rear (~1.7 m / 67 in) | 67 in | 55 in | 2.35 × 55 × 55 ≈ **7,100 µF** |

Because every run exceeds 3 feet, this capacitance cannot sit only at the ESC. It must be spread in stations every 8 to 12 inches along each boom, roughly 5 to 7 stations per boom, with the final bank within 4 inches of the ESC.

Practical shape of the solution:

- Low ESR electrolytics (Rubycon ZLH, Panasonic FR, Nichicon PW class or equivalent). 63 V is the minimum rating over the 50.4 V full-charge pack and leaves only 20 percent headroom against exactly the spikes the caps exist to absorb. 80 V is the comfortable rating where can size allows.
- Roughly 6 to 7 × 1,000 µF cans per boom, about 25,000 µF fleet total. Estimated added mass ~250 to 400 g across the aircraft including boards, leads, and mounting, before any potting or strain relief.
- Electrolytic cans mounted along booms sit in vibration and prop wash for their whole life. Every solder joint on the HV bus is a new failure point, and the team has flagged that many purchased components are last in stock, so rework tolerance is low.
- The full bank adds to QS8-S connection inrush. The anti-spark resistor handles it, but precharge lengthens.
- Warranty on the ESCs is void per KDE's own wording.

What Option A buys: the ESC heatsink sits in prop wash (well above the 5 mph requirement in hover), the motor leads stay factory length, and the mount design continues unchanged.

## Option B: ESC in the fuselage, extend the motor (AC) leads

This is KDE's recommended configuration.

- DC path: FCHUB pads to each ESC inside the CG bay, well under 12 inches. No added capacitance, no warranty issue.
- AC path: each motor's three 14 AWG phase leads extend ~1.4 to 1.7 m along the boom. Extending motor leads is the configuration KDE explicitly recommends, and phase lead inductance is negligible against the motor winding inductance at these lengths. Twist or braid each motor's three phase wires along the boom to cut radiated switching noise, and keep separation from the CAN pair being pulled for the future MAD arrangement.
- The 22 AWG PWM control lead stays inside the fuselage at the flight controller. This eliminates the PWM signal integrity risk over boom-length runs beside HV cabling that was flagged on the July 20 call. The spare conductor should still be pulled with the boom CAN cable so the MAD DroneCAN arrangement can return without re-pulling booms.
- Mass: three 14 AWG phase wires weigh slightly more per meter than two 13 AWG DC wires (~25 to 30 g more per boom, ~110 g fleet), but no capacitor banks and no boom PWM leads. Net roughly 200 g lighter than Option A.
- Thermal is the one real cost. Rough loss estimate at 30 to 40 A hover is on the order of 10 to 20 W per ESC, 40 to 80 W in the bay plus regulator losses. The 55 A rating assumes 5 mph airflow, and a closed bay has none. Mitigations, in order of preference: bolt the ESC heatsinks to an aluminum plate (the KDE FAQ allows a metal plate as the heatsink), add intake and exhaust vents to the bay, add a small fan if needed. The two I2C bay temperature sensors already planned (July 10 call) answer the fan question with flight data. Stork's fuselage-mounted 4-in-1 is the in-house precedent.

## Option C: ESC at the boom root, just outside the skin

A compromise worth checking against the actual structure. The FCHUB sits at the boom root convergence, so an ESC mounted on each boom root outside the skin may keep the DC run under 12 inches while placing the heatsink in external air and prop wash.

- No capacitors, short DC, short control leads, free airflow.
- AC extension shortens to roughly 1.2 to 1.5 m.
- Depends entirely on whether ~12 inches of DC reach is physically achievable from the FCHUB pads to a mountable root position, and on clearance with the wing to fuselage transition pieces. If the DC run creeps past 12 inches this collapses into a small version of Option A.

## Comparison

| | A: at motor | B: in fuselage | C: boom root |
|---|---|---|---|
| KDE position | Own risk, warranty void | Recommended | Compliant if under 12 in |
| Capacitors | ~5,600 to 7,100 µF per ESC, distributed along booms | None | None |
| Added mass | ~250 to 400 g | ~110 g wire | ~80 g wire |
| PWM lead exposure | Full boom beside HV | None | Short |
| Cooling | Prop wash, excellent | Bay, needs a plate/vents/fan | External air, good |
| Build complexity | High (cap stations, HV joints) | Low | Low, geometry permitting |
| Failure points added | Many | Few | Few |
| Mount design impact | None | Rework | Rework |

## Recommendation

**Option B, ESCs in the fuselage with extended, twisted motor leads.** It is the vendor's recommended configuration, removes the capacitor engineering entirely, deletes the PWM over boom risk, is the lightest, and is reversible when the MAD DroneCAN set arrives (booms carry phase wires plus the pre-pulled CAN and spare conductor either way). Fall back to Option C if bay thermal proves unworkable and the boom root geometry supports a sub 12 inch DC run. Option A should only survive if a structural or mount constraint rules both others out, in which case the capacitor specification above applies as written.

Condition attached to the recommendation: bench-run one KDE motor and ESC on the testbed at hover-representative load and log current and ESC temperature. That single test validates both the thrust ceiling question from the July 20 call and the bay heat budget this recommendation depends on.

## Impact if Option B is adopted

- SPH-E-001 H3 changes for PT1: the 12 AWG DC extensions and XT90S boom junction splices are not used with the KDE set. Boom runs carry twisted 14 AWG phase extensions instead. The purchased XT90 stock stays for the MAD return and the phase extension joints can use KDE 3.5/4.0 mm bullets or soldered splices.
- The tail HV spur (H5) breakout point moves: with no ESC lead extension joint at the tail boom root, the 18 AWG tail feed taps the FCHUB side directly. Small harness edit, no new parts.
- CG bay layout gains four ESC positions plus a heatsink plate. Compartment is ~180 × 250 × 90 mm, ESC body is 35 × 72 mm, so four fit with the plate setting the packaging.
- Parts list additions: 14 AWG silicone wire lot, bullet connectors or splice sleeves for twelve phase joints, optional bay fan. Capacitor line items not needed.

## Sources

- KDE-UAS55HVC product page: https://www.kdedirect.com/collections/uas-multi-rotor-electronics/products/kde-uas55hvc
- KDE Direct, The Dangers of Long DC Wires: https://www.kdedirect.com/blogs/news/dc-wire-length
- KDE Direct Support Center FAQ, Electronics: https://www.kdedirect.com/pages/support-center-faq-electronics
- SPH-E-001 §4.14 harness table (run lengths), §7 thermal notes
