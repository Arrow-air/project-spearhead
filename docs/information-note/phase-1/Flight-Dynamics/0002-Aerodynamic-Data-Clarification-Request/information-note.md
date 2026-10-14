---
title: "Aerodynamic Data Clarification Request"
sidebar_label: "Information Note"
sidebar_position: 1
---

# Aerodynamic Data Clarification Request

## Status

"Open — awaiting Aero team response"

"Revision History: Rev A, 14 October 2026, initial release (flight dynamics, ZeynepB5793)"

"Replacement Log: None"

"Reference: [Preliminary Design Aerodynamic Database](../../Aerodynamics/0002-Preliminary-Design-Aerodynamic-Database/information-note.md), [Preliminary Fixed-Wing Stability Analysis](../0001-Preliminary-Fixed-Wing-Stability-Analysis/information-note.md)"

## Purpose

The October aerodynamic database (`tail_update_2026_10_07` base table and the 11 October 103-case control delivery) is being used for flight dynamics work. Before the flight mechanics model is corrected, a few reference, geometry and data points need confirmation from the Aero team. CG location and pusher motor selection will be decided afterwards, jointly, using the confirmed data.

The questions are limited to what can be answered within one working day from existing CAD, setup files, exports and results. The only new CFD request is optional (AQ-12). Please answer in the [Response Log](#response-log) or reply per question ID. If a question rests on a wrong reading of the data, say so; that is also an answer.

Priorities:

- **P0** — needed before the model can be corrected.
- **P1** — useful for the same correction pass; answer if time allows.

## Conventions Taken From the Deliveries

Please confirm or correct each line.

| Item | Our reading |
|---|---|
| Body axes | x forward, y right, z down; right-hand moments |
| Moment reference of the base and control tables | (0, 0, 0) m, the raw CFD origin (confirmed from 6 October) |
| Physical location of the CFD origin | Wing quarter-chord on the chord line at the root; inferred from the preliminary STEP files, not confirmed |
| Preliminary CAD frame | x aft, z up, same origin as the CFD (inferred) |
| Reference values | S = 1.805 m², b = 3.95 m, c = 0.457 m, V = 25 m/s, rho = 1.021 kg/m³ |
| Base export header | Records final moment center (−0.1, 0, +0.15) m, which conflicts with the confirmed origin |
| Corrections already in the base table | 1.3× wind-axis drag multiplier; Fy, Mx, Mz zeroed at beta = 0 |
| Surface mapping proposed by the Aero team | Left/right aileron = +da/−da; left/right ruddervator = de+dr / de−dr |
| Proposed surface travel | ±30° for all surfaces (11 October) |

## Values Derived From the October Data

These were calculated from the delivered tables. Please confirm them from your own post-processing, or point out where our reading differs.

| Quantity | Our value |
|---|---|
| Stick-fixed neutral point, alpha 0–8°, 25 m/s | x ≈ −0.05 m (z reference 0) to −0.07 m (z reference +0.15 m); about 36–40 % MAC if the origin is the wing quarter-chord |
| Pitching moment at the origin, alpha = beta = 0, 25 m/s | My = −37.54 N·m (Cm = −0.143); Cm ≈ −0.07 at CL = 0 |
| Lift-curve slope, alpha 0–8° | CL_alpha ≈ 5.07 /rad |
| Maximum lift | CLmax ≈ 1.386 at alpha = 16° |
| Directional and lateral static stability (moments near the origin) | Cn_beta ≈ +0.03 to +0.04 /rad; Cl_beta ≈ −0.002 to −0.020 /rad depending on alpha and reference height |
| Neutral-row My difference between base and control tables at alpha = beta = 0, 25 m/s, after restoring the base to the origin | 0.972 N·m |
| Elevator coverage in the 11 October delivery | −15° to +15° |

## Questions

### A. Geometry and Reference Frame

#### AQ-01 — Physical location of the CFD origin (P0)

- **Question:** Where exactly is (0, 0, 0) on the CFD geometry? Is it the wing quarter-chord on the root chord line, the same origin as the preliminary STEP files (`main_wing.step`, `tail_wing.step`, `fuselage.step`, `longitudinal_booms.step`)? Which axis directions does the CFD geometry use?
- **Deliverable:** one sentence, plus a screenshot of the CFD geometry with the origin marked.

#### AQ-02 — Moment reference of each delivered table (P0)

- **Question:** From the Nondimit export settings, which moment reference does (a) the October base export and (b) the 11 October control table actually use? The base header says (−0.1, 0, +0.15) m while the confirmed origin is (0, 0, 0). After restoring the base numerically to (0, 0, 0), the neutral rows of the two tables still differ by 0.972 N·m in My at alpha = beta = 0. Is that residual expected (surrogate fit, geometry revision, correction difference)?

#### AQ-03 — October tail geometry (P0)

- **Question:** Please supply the October tail as meshed: root and tip chord, planform area, geometric and projected span, dihedral angle, incidence, airfoil (the design note says NACA 0015, the sizing configuration says NACA 0018), ruddervator chord fraction and hinge line, and the x, z position of the tail root and apex.
- **Deliverable:** a table, plus the STEP of the meshed geometry revision if it is at hand (earlier case paths mention `v14`).

#### AQ-04 — Wing and tail incidence and sign convention (P0)

- **Question:** Confirm the wing incidence (+1°) and tail incidence (+2°) used in the CFD geometry, and the sign convention (positive = leading edge up relative to the fuselage datum?).

#### AQ-05 — Components in the mesh and force integration (P0)

- **Question:** Which components were meshed and integrated: fuselage, booms, lift motors, stopped lift propellers, landing gear, pusher nacelle? Are the actuator-disk forces separated from the airframe loads in the delivered tables? What does the 1.3 drag multiplier represent, and is it applied identically in the base and control tables?
- **Deliverable:** list of integrated surfaces per dataset.

#### AQ-06 — Booms and lift motor pods in the CFD geometry (P1)

- **Question:** Are the booms and lift motor pods in the CFD geometry? If yes, at which x, y, z positions relative to the origin?

### B. Longitudinal Stability

#### AQ-07 — Neutral point (P0)

- **Question:** Do you agree that the October database puts the stick-fixed neutral point at x ≈ −0.05 to −0.07 m (about 36–40 % MAC)?
- **Deliverable:** the neutral point from your own post-processing of the existing data.

#### AQ-08 — Tail effectiveness and fuselage contribution (P1)

- **Question:** A vortex-lattice check of the wing and tail alone (no fuselage, booms or nacelles) puts the neutral point at about 47–49 % MAC, against 36–40 % from the CFD. From the existing flow fields, is a fuselage, boom and nacelle contribution of about 0.04–0.06 m plausible, and is the tail in the wing wake at cruise alpha? A qualitative answer is enough.

### C. Control Surfaces

#### AQ-09 — Ruddervator sign convention and hinge angles (P0)

- **Question:** The delivery uses virtual elevator and rudder columns. For the existing cases, what physical hinge angle did each ruddervator have, and what is the sign convention? For this tail (apex above the booms) we expect a nose-right yaw moment from the left trailing edge down and the right trailing edge up, the opposite of an upright V. Is that correct?

#### AQ-10 — Wing surfaces: aileron or flaperon (P1)

- **Question:** The CFD models plain ailerons; the electrical master document calls the wing surfaces flaperons. Which is the design? Please give span station, chord fraction and hinge line of the modeled surfaces.

#### AQ-11 — Hinge moments from existing cases (P1)

- **Question:** If each control surface is a separate patch in the existing meshes, please export the hinge moments of the ruddervators and ailerons for the existing ±15° (elevator, rudder) and ±30° (aileron) cases at 25 m/s. Skip if it needs a re-run.
- **Why:** servo check against Kingmax CLS3015S, 35 kg·cm at 8.4 V.

#### AQ-12 — Minimal elevator extension (P1, optional new CFD)

- **Question:** Only if it fits in the day: elevator at −20° and −30° using the existing deflected-elevator setup.

| Speed (m/s) | Alpha (°) | Beta (°) | Elevator (°) | Cases |
|---|---|---|---|---:|
| 22 | 0, 2, 4 | 0 | −20, −30 | 6 |

Moments about (0, 0, 0); dimensional loads Fx, Fy, Fz, Mx, My, Mz; lift rotors off; same geometry revision as the 11 October delivery.

### D. Validity Range of the Existing Data

#### AQ-13 — Speed and stall (P1)

- **Question:** The dense base table is at 25 m/s only. In your judgement, are the coefficients usable at 15–20 m/s (Re ≈ 0.4–0.6 million)? How confident is the stall (CLmax ≈ 1.386 at 16°) in steady RANS?

#### AQ-14 — Large sideslip and controls in sideslip (P1)

- **Question:** Is the base table reliable at beta = 30–90° and low speed, given that it comes from a 25 m/s surrogate? Can the beta = 0 control increments be used at small sideslip (±5–10°)?

### E. Delivery Documentation

#### AQ-15 — Convergence evidence (P1)

- **Question:** If the solver logs from the October runs are kept, please share per-case residual and force-history plots, plus solver, turbulence model, mesh size, y+ and iteration limit. Earlier metadata says `status=ok`, `convergence=not_assessed`, `max_iterations=400`.

#### AQ-16 — Manifest for each delivery (P1)

- **Question:** For each delivered file, a short manifest: geometry revision, reference values, moment center, density, velocity, corrections applied (drag multiplier, beta = 0 zeroing), column definitions and signs. The 11 October control table has no embedded metadata and no coefficient columns.

#### AQ-17 — Design density (P1)

- **Question:** Confirm that rho = 1.021 kg/m³ represents the Ankara test field (about 1125 m) on a hot day (about 30 °C). On a 15 °C day it is closer to 1.08. Is 1.021 intended as the conservative design case for all analyses?

## Response Log

| ID | Priority | Answer | Answered by | Date |
|---|---|---|---|---|
| AQ-01 | P0 | | | |
| AQ-02 | P0 | Confirmed: `database.csv` Nondimit export center was (−0.1, 0, +0.15); deflection tables are at (0, 0, 0). The 0.972 N·m residual remains open; the drag multiplier is applied before the shift, so it does not explain it. | Aero team | 14 Oct 2026 |
| AQ-03 | P0 | | | |
| AQ-04 | P0 | | | |
| AQ-05 | P0 | | | |
| AQ-06 | P1 | | | |
| AQ-07 | P0 | | | |
| AQ-08 | P1 | | | |
| AQ-09 | P0 | Sign conventions confirmed | Aero team | 14 Oct 2026 |
| AQ-10 | P1 | | | |
| AQ-11 | P1 | | | |
| AQ-12 | P1 | | | |
| AQ-13 | P1 | | | |
| AQ-14 | P1 | | | |
| AQ-15 | P1 | | | |
| AQ-16 | P1 | | | |
| AQ-17 | P1 | | | |
