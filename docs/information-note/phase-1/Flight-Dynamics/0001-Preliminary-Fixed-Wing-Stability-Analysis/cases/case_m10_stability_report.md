# Case M10 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_m10.csv`
- Moment/reference/CG: `(-0.100, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `1.7432 deg`
- theta: `1.7432 deg`
- elevator: `-1.8846 deg`
- throttle: `0.1608`
- residual norm: `5.859e-14`

## Static Stability Derivatives

- Cm_alpha: `+0.446907 /rad` (unstable longitudinal static sign)
- Cn_beta: `+0.028648 /rad` (stable directional static sign)
- Cl_beta: `-0.010723 /rad` (stable lateral static sign)
- CL_alpha: `+4.858682 /rad`
- CY_beta: `-0.227790 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| phugoid | `-0.32842+0.893381j; -0.32842-0.893381j` | 0.9518 | 0.3450 | 3.0449 | stable |
| short_period | `-9.31431+0j; +1.1421+0j` | 1.1421 | -1.0000 | 0.8756 | unstable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.2281+0j` | 30.2281 | 1.0000 | 0.0331 | stable |
| spiral | `+0.000287463+0j` | 0.0003 | -1.0000 | 3478.7093 | unstable |

## Summary

- Static stability: longitudinal `unstable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `unstable`.
- Unstable modes: short_period, spiral.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case M10 improves the longitudinal trend relative to nominal but remains unsuitable as a stability baseline. The pitch static derivative is still unstable, the longitudinal dynamic response is unstable, and the spiral mode is weakly unstable in this linearization.
