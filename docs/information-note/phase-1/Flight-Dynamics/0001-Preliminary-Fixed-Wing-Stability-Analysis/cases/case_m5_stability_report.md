# Case M5 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_m5.csv`
- Moment/reference/CG: `(-0.050, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `2.0061 deg`
- theta: `2.0061 deg`
- elevator: `-5.5737 deg`
- throttle: `0.1665`
- residual norm: `4.705e-14`

## Static Stability Derivatives

- Cm_alpha: `-0.291354 /rad` (stable longitudinal static sign)
- Cn_beta: `+0.028665 /rad` (stable directional static sign)
- Cl_beta: `-0.011511 /rad` (stable lateral static sign)
- CL_alpha: `+4.932043 /rad`
- CY_beta: `-0.226353 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| phugoid | `-0.024601+0.474923j; -0.024601-0.474923j` | 0.4756 | 0.0517 | 40.6487 | stable |
| short_period | `-4.45504+5.15473j; -4.45504-5.15473j` | 6.8131 | 0.6539 | 0.2245 | stable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.2267+0j` | 30.2267 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00364264+0j` | 0.0036 | 1.0000 | 274.5265 | stable |

## Summary

- Static stability: longitudinal `stable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `stable`.
- Weakly damped modes: phugoid.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case M5 is the recommended near-term development baseline from this sweep. It preserves stable identified longitudinal behavior with more margin than M5p5 while staying near the transition region that should be refined in future CG studies.
