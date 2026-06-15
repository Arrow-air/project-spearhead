# Case 0 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_0.csv`
- Moment/reference/CG: `(-0.000, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `2.2937 deg`
- theta: `2.2937 deg`
- elevator: `-9.7286 deg`
- throttle: `0.1719`
- residual norm: `1.415e-14`

## Static Stability Derivatives

- Cm_alpha: `-1.037054 /rad` (stable longitudinal static sign)
- Cn_beta: `+0.034377 /rad` (stable directional static sign)
- Cl_beta: `-0.013984 /rad` (stable lateral static sign)
- CL_alpha: `+5.001922 /rad`
- CY_beta: `-0.228001 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| phugoid | `-0.0253188+0.540732j; -0.0253188-0.540732j` | 0.5413 | 0.0468 | 39.4963 | stable |
| short_period | `-4.45612+7.50931j; -4.45612-7.50931j` | 8.7319 | 0.5103 | 0.2244 | stable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.2252+0j` | 30.2252 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00285347+0j` | 0.0029 | 1.0000 | 350.4503 | stable |

## Summary

- Static stability: longitudinal `stable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `stable`.
- Weakly damped modes: phugoid.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case 0 is the most stable sampled case by longitudinal static and dynamic metrics, but it requires the largest elevator trim deflection. It should be treated as a forward stability bound until control authority, mass properties, and the full control model are validated.
