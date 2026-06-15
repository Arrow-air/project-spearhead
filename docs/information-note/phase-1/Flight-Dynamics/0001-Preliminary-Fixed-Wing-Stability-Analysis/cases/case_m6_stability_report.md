# Case M6 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_m6.csv`
- Moment/reference/CG: `(-0.060, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `1.9516 deg`
- theta: `1.9516 deg`
- elevator: `-4.8077 deg`
- throttle: `0.1653`
- residual norm: `2.371e-14`

## Static Stability Derivatives

- Cm_alpha: `-0.138827 /rad` (stable longitudinal static sign)
- Cn_beta: `+0.028648 /rad` (stable directional static sign)
- Cl_beta: `-0.011321 /rad` (stable lateral static sign)
- CL_alpha: `+4.916448 /rad`
- CY_beta: `-0.226595 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| short_period | `-5.50085+0j; -3.20677+0j; -0.240531+0j; +0.116444+0j` | 0.1164 | -1.0000 | 8.5878 | unstable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.227+0j` | 30.2270 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00266744+0j` | 0.0027 | 1.0000 | 374.8912 | stable |

## Summary

- Static stability: longitudinal `stable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `unstable`.
- Unstable modes: short_period.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case M6 is statically stable in pitch but still dynamically unstable. It sits close to the dynamic stability transition and is useful for bounding the aft side of the stable region.
