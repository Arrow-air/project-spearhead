# Propulsion Model

Engine and propeller performance model for the Spearhead UAV. Computes power, thrust, fuel consumption, endurance, and transition (acceleration to cruise) for a 2-stroke gasoline engine with fixed-pitch propeller.

## Quick Start

```bash
python propulsion_model.py
```

Generates an 8-panel figure (`spearhead_engine_model.png`) and prints a console summary with cruise RPM, fuel flow, endurance, rate of climb, and prop comparison table.

## Configuration

All inputs are in the `DEFAULT_CONFIG` dict at the top of `propulsion_model.py`. Pass overrides to `make_plots()` or edit the dict directly.

## Outputs

### Figure (8 subplots)
| Plot | Description |
|------|-------------|
| Power & Torque vs RPM | Shaft power and torque curves from VE model |
| Static Thrust vs RPM | Multiple prop candidates, anchored to manufacturer data |
| Fuel Consumption vs RPM | Fuel flow from BSFC model |
| Endurance vs Tank Size | Tank sizing for target endurance at cruise |
| Force Balance vs Speed | Thrust and drag curves, net force |
| Airspeed vs Time | Transition simulation (acceleration to cruise) |
| Acceleration vs Speed | Net acceleration during transition |
| Distance vs Time | Ground distance covered during transition |

### Console Summary
- Cruise RPM (where thrust = drag)
- Fuel flow at cruise
- Excess power and max rate of climb
- Battery charging time from excess shaft power
- Prop comparison table across all candidates

## Models and Equations

### Engine Power
```
P = P_rated * (VE(RPM) / VE(RPM_rated)) * (RPM / RPM_rated)
```
- VE shape: quadratic + cubic around peak, asymmetric falloff (faster drop above peak)
- `VE = VE_peak * (1 - 0.55*(x-1)^2 - 0.07*(x-1)^3)` where `x = RPM/RPM_ve_peak`
- Anchored to manufacturer rated HP at rated RPM

### Torque
```
Q = P / omega = (P_kW * 1000) / (RPM * 2*pi/60)
```

### Static Thrust
```
T_static = T_anchor * prop_scale * (rho/rho_anchor) * (RPM/RPM_anchor)^(4/3)
```
- RPM^(4/3) exponent from combined P ~ RPM^2 (2-stroke) and momentum theory T ~ P^(2/3)
- Density ratio referenced to anchor test altitude (not sea level)
- Prop scaling: `(D/D_ref)^(4/3) * (p/p_ref)^(1/3)`

### Thrust in Flight
```
T = T_static * max(0, 1 - (J/J_zero)^2)
```
- Advance ratio: `J = V / (n * D)` where `n = RPM/60`
- Zero-thrust advance ratio: `J_zero = eta_pitch * (p/D)`
- `eta_pitch` = 1.10-1.30 for fixed-pitch props (thrust extends beyond geometric pitch speed)

### Propeller Efficiency
```
eta = eta_max * 4 * x * (1-x)   where x = J/J_zero
```
- Parabolic model peaking at J = 0.5 * J_zero with value eta_max
- eta = 0 at static (J=0) and at zero-thrust (J=J_zero)

### Fuel Consumption
```
fuel_flow [ml/hr] = P_shaft [kW] * BSFC [g/kWh] / 0.74 [kg/L gasoline]
```

### Drag
```
D = D_cruise * (V / V_cruise)^2
```
- Simple quadratic scaling from a single calibration point (from sizing tool)

### Transition Simulation
- 1D horizontal Euler integration, dt = 0.02 s
- `F_net = T(V, RPM) - D(V)`, `a = F_net / mass`, `V += a * dt`

## Assumptions

### Engine
- 2-stroke piston-port timing: VE curve is asymmetric with faster falloff above peak RPM
- Power scales linearly with VE and RPM: `P ~ VE * RPM` (constant displacement)
- BSFC is constant across the operating range (520 g/kWh default)
- No altitude derating on power — only thrust is density-corrected

### Propeller
- Fixed pitch — no variable-pitch or constant-speed governor
- Ct(J) follows a quadratic decay: `1 - (J/J_zero)^2`
- Pitch scaling `(p/p_ref)^(1/3)` is empirical, not derived from BEM theory
- Diameter scaling `(D/D_ref)^(4/3)` assumes engine absorbs prop load at same RPM
- eta_pitch > 1 accounts for real props producing thrust beyond geometric pitch speed
- Propeller efficiency model is a symmetric parabola — real props have an asymmetric eta(J) curve

### Thrust
- Static thrust anchored to a single manufacturer data point (9 kg @ 7900 RPM, 20x10 prop)
- RPM^(4/3) scaling assumes P ~ RPM^2, which is approximate for 2-strokes above torque peak
- Density correction is linear: `T ~ rho` (momentum theory gives `T ~ rho^(1/3)` at constant power, but combined engine+prop gives ~linear)
- No compressibility, no tip Mach effects

### Drag
- Quadratic V^2 scaling from a single point — does not capture the induced drag polar shape
- No density correction on drag (drag_cruise_N is assumed at mission altitude)

### Transition
- Level flight only — no ground roll friction, no climb, no wind
- Constant RPM throughout transition (no throttle modulation)
- Euler integration (first-order, dt = 0.02 s)

## Dependencies

```
numpy
matplotlib
```

## Calibration Data

| Parameter | Value | Source |
|-----------|-------|--------|
| Rated HP | 4.1 HP @ 9000 RPM | RCGF Stinger 35cc RE datasheet |
| VE peak RPM | 7500 RPM | Typical for piston-port 2-stroke |
| BSFC | 520 g/kWh | Typical for small 2-stroke (not measured) |
| Anchor thrust | 9 kg @ 7900 RPM | Manufacturer test, 20x10 prop, ~100 m altitude |
| eta_pitch | 1.20 | Typical for fixed-pitch RC props (1.10-1.30) |
| eta_prop_max | 0.70 | Typical for fixed-pitch (0.60-0.75) |
| D_cruise | 30 N @ 25 m/s | From initial sizing tool |

## Validation

### Thermal Efficiency
The model's BSFC of 520 g/kWh yields a thermal efficiency of 16.0%, calculated as:

```
eta_th = 3600 / (BSFC [g/kWh] * LHV [kJ/g])
       = 3600 / (520 * 43.4)
       = 16.0%
```

This is consistent with published data for small piston-port 2-stroke engines:

| Engine class | BSFC (g/kWh) | Thermal eff. |
|-------------|-------------|-------------|
| Small RC 2-stroke (piston-port) | 450 - 600 | 14 - 18% |
| Tuned 2-stroke (DI / stratified) | 320 - 450 | 18 - 26% |
| Small 4-stroke (Saito, OS) | 290 - 370 | 22 - 28% |
| Automotive gasoline | 240 - 280 | 30 - 35% |

The RCGF Stinger 35cc RE is a piston-port 2-stroke with no fuel injection or tuned exhaust, placing it squarely in the 14-18% range. No dyno data is available for this specific engine; the BSFC of 520 g/kWh is an estimate based on the engine class, not a measurement.

### Thrust Anchor
Static thrust is anchored to a single manufacturer test point (9 kg at 7900 RPM with a 20x10 prop at ~100 m altitude). The density ratio correction uses the anchor altitude as reference, so the anchor point reproduces exactly at its test conditions. Thrust at other RPMs and altitudes is extrapolated via the RPM^(4/3) scaling law.

### Cruise RPM
The cruise RPM solver finds the operating point where thrust equals drag at 25 m/s. For the 20x10 prop at 2000 m altitude, this yields ~7000 RPM. No flight test data is available for comparison.

### Limitations
- BSFC is assumed constant across all RPMs. In reality it varies by 20-30% between peak torque and maximum RPM.
- The thrust model has no compressibility or tip Mach correction — valid for tip speeds below ~200 m/s (Mach 0.6).
- The Ct(J) quadratic model is a simplification of the real thrust coefficient curve, which is typically more linear in the mid-J range.
- Propeller efficiency uses a symmetric parabolic model; real eta(J) curves are asymmetric with a sharper cutoff near J_zero.
