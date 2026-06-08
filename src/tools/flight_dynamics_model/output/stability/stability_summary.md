# Preliminary Stability Summary

> **PRELIMINARY / NOT VALIDATED: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.**

## Trim

- Success: `True`
- Alpha: `-0.0347 deg`
- Theta: `-0.0347 deg`
- Elevator: `20.7092 deg`
- Throttle: `0.13543`
- Residuals: `[ 4.77396e-17 -1.13687e-15 -1.66533e-15 -1.08420e-19]`

## Static Derivatives

- `CL_alpha = +4.448455 /rad`
- `Cm_alpha = +5.195033 /rad`
- `Cl_beta = +0.002865 /rad`
- `Cn_beta = -0.003424 /rad`
- Windowed diagnostic slopes over `+/-5.0 deg`: `CL_alpha = +4.368155`, `Cm_alpha = +5.143861`, `Cl_beta = -0.003964`, `Cn_beta = -0.001067` `/rad`
- Inferred neutral-point station relative to the assumed CG, body x-forward: `+0.5337 m`

## Longitudinal Modes

- `non-oscillatory          lambda=-22.46121+0.00000j wn=22.4612 rad/s zeta=1.0000`
- `non-oscillatory          lambda=+13.00573+0.00000j wn=13.0057 rad/s zeta=-1.0000`
- `phugoid                  lambda=-0.03510+0.65895j wn=0.6599 rad/s zeta=0.0532`
- `phugoid                  lambda=-0.03510-0.65895j wn=0.6599 rad/s zeta=0.0532`

## Lateral-Directional Modes

- `heading integrator       lambda=+0.00000+0.00000j wn=0.0000 rad/s zeta=n/a`
- `roll subsidence          lambda=-33.96210+0.00000j wn=33.9621 rad/s zeta=1.0000`
- `yaw/directional-related  lambda=-9.83677+0.00000j wn=9.8368 rad/s zeta=1.0000`
- `spiral                   lambda=+0.11569+0.00000j wn=0.1157 rad/s zeta=-1.0000`
- `yaw/directional-related  lambda=-0.16401+0.00000j wn=0.1640 rad/s zeta=1.0000`

## Limitations

- The inertia matrix is a placeholder.
- ADB v1.1 contains alpha-beta data only.
- Control increments and rate damping derivatives are provisional placeholders.
- Propulsion is a placeholder pusher thrust model.
