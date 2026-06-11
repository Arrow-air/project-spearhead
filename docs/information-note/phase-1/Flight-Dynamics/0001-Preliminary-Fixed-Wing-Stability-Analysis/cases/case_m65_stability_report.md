# Case M65 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_m65.csv`
- Moment/reference/CG: `(-0.065, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `1.9248 deg`
- theta: `1.9248 deg`
- elevator: `-4.4307 deg`
- throttle: `0.1647`
- residual norm: `2.276e-14`

## Static Stability Derivatives

- Cm_alpha: `-0.065699 /rad` (stable longitudinal static sign)
- Cn_beta: `+0.028648 /rad` (stable directional static sign)
- Cl_beta: `-0.011244 /rad` (stable lateral static sign)
- CL_alpha: `+4.908748 /rad`
- CY_beta: `-0.226749 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| short_period | `-6.39647+0j; -2.24543+0j; -0.442251+0j; +0.252802+0j` | 0.2528 | -1.0000 | 3.9557 | unstable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.2271+0j` | 30.2271 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00228627+0j` | 0.0023 | 1.0000 | 437.3945 | stable |

## Summary

- Static stability: longitudinal `stable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `unstable`.
- Unstable modes: short_period.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case M65 is the first sampled case with the stable sign for pitch static stability, but it is not dynamically acceptable. It helps bracket the transition region between static stability and fully stable identified longitudinal dynamics.
