# ArduPilot QuadPlane Reference for Project Spearhead

*Extracted and filtered from the official ArduPilot Plane documentation. Tailsitter, tilt-rotor, and other non-applicable sections omitted. Focused on a Quad-X pusher IC engine configuration at ~25 kg MTOW.*

---

## 1. Enabling QuadPlane

QuadPlane is a mode within the standard **Plane** firmware — there is no separate firmware. Set `Q_ENABLE = 1` and reboot. All QuadPlane parameters use the `Q_` prefix. After enabling, refresh the parameter list to see the full set.

---

## 2. Frame Setup

### Frame Class & Type

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `Q_FRAME_CLASS` | 1 | Quad (4 VTOL motors) |
| `Q_FRAME_TYPE` | 1 | X frame layout |

For Spearhead's quad-X, use `Q_FRAME_CLASS = 1`, `Q_FRAME_TYPE = 1`.

**H vs X mixing:** If the airframe is perfectly rigid the difference is negligible. However, quadplane structures are often less rigid than traditional multirotors and torsional effects must be considered. When yaw input results in most weight being supported by a diagonal motor pair, the fuselage or wing may twist and distort thrust vectors. Choose the mixing scheme (X or H) that ensures any induced twist *complements* the desired yaw. For motors mounted on arms extending fore/aft aligned with the wing chord, X mixing is usually correct.

### Motor Ordering (Quad-X Default)

| Output | Motor | Position | Rotation |
|--------|-------|----------|----------|
| 5 | Motor 1 | Front right | Counter-clockwise |
| 6 | Motor 2 | Rear left | Counter-clockwise |
| 7 | Motor 3 | Front left | Clockwise |
| 8 | Motor 4 | Rear right | Clockwise |

Rule of thumb: "motors turn in towards the fuselage" for X frame.

Outputs 1–4 are reserved for conventional plane control surfaces (aileron, elevator, rudder, throttle) by default when `Q_ENABLE` is set.

### Remapping Outputs

If you need motors on different channels, set `SERVOn_FUNCTION` values: Motor 1 = 33, Motor 2 = 34, Motor 3 = 35, Motor 4 = 36. Only do this for non-standard ordering; the defaults auto-configure on boot.

---

## 3. Building Guidance

### Key Design Principles

- **Frame strength:** Must carry VTOL motors, power system, and payload without flex. Minimum wing/frame/mount twist is critical so motors provide thrust vertically at all times.
- **Motor clearance:** Complete clearance above and below the full disk area of each VTOL motor for full aerodynamic thrust.
- **Robust mounts:** The VTOL motor mounting system must handle the vibration and loads of hover flight at full MTOW.
- **Drag minimization:** Minimize aerodynamic drag from VTOL motors, arms, and mounting hardware for efficient cruise.
- **Motor alignment is critical:** A few degrees of misalignment can eliminate yaw authority in one or both directions. Check on the bench by blocking the vehicle level and measuring prop-tip-to-table distances in fore-aft and side-to-side orientations. The arc-sine of tip-end differences divided by prop diameter gives the tilt angle. Even 1–2° matters for yaw.
- **Intentional yaw tilt:** You can *increase* yaw authority by purposely tilting one or both pairs of adjacent rotating motors 1–2° in the direction of their torque (outward for X frame).
- **Use eCalc** to choose motors, ESCs, batteries, and propellers. Motor power-to-weight ratios vary enormously.

### Wing Downforce in Hover

The VTOL motors must provide thrust not just for total airframe weight, but also for the additional downforce on the wings when hovering in a horizontal stance. This is a significant consideration at 25 kg MTOW with a 3.95 m span wing.

### IC Engine for Cruise — Range Advantage

A QuadPlane with VTOL takeoff can actually have *greater* range than the same airframe with conventional launch, because:

1. No need for high-thrust takeoff capability on the forward motor — allows optimizing the pusher prop/engine for cruise efficiency (larger prop, geared motor, etc.)
2. No stall-speed acceleration requirement at launch — can carry more fuel/battery
3. IC engine cruise provides dramatically better endurance than electric

To maximize this advantage, do very rapid VTOL takeoffs and landings to minimize battery consumption in hover.

---

## 4. Critical Parameters

### Core QuadPlane Parameters

| Parameter | Purpose | Spearhead Notes |
|-----------|---------|-----------------|
| `Q_ENABLE` | Enable QuadPlane (1) | Must set and reboot |
| `Q_FRAME_CLASS` | Motor arrangement | 1 (Quad) |
| `Q_FRAME_TYPE` | Motor layout | 1 (X) |
| `Q_M_PWM_MIN` / `Q_M_PWM_MAX` | PWM range for VTOL ESCs | Set to match your ESC protocol |
| `Q_M_SPIN_ARM` | Motor output when armed in quad mode | Important for idle behavior |
| `Q_M_THST_HOVER` | Throttle % for hover at mid-stick | Calibrate so QSTABILIZE hovers at mid-stick ±6% |
| `Q_M_HOVER_LEARN` | Auto-learn hover throttle | Enable in QLOITER/QHOVER |
| `Q_A_RAT_RLL_P` | Roll rate P gain | Default 0.25 — QuadPlanes often need significantly higher |
| `Q_A_RAT_PIT_P` | Pitch rate P gain | Default 0.25 — often needs to be higher |
| `Q_A_ANGLE_MAX` | Max lean angle in VTOL modes (cdeg) | Controls pitch/roll limits |
| `Q_TRIM_PITCH` | VTOL pitch trim offset | Corrects backward drift from fixed-wing level trim |
| `SCHED_LOOP_RATE` | Scheduler loop rate | Default 300 Hz for QuadPlane — do not raise for 25 kg vehicle |
| `ARMING_RUDDER` | Rudder arming | Set to 2 for rudder disarm |

### Transition Parameters

| Parameter | Purpose | Spearhead Notes |
|-----------|---------|-----------------|
| `AIRSPEED_MIN` | Minimum airspeed / transition complete speed | VTOL motors cut after reaching this |
| `Q_TRANSITION_MS` | Time to ramp down VTOL motors after reaching AIRSPEED_MIN | Default 5000 ms (5 sec) |
| `Q_TRANS_FAIL` | Transition timeout (sec) | 0 = disabled. If set, transition aborts if AIRSPEED_MIN not reached |
| `Q_TRANS_FAIL_ACT` | Action on transition failure | What happens when Q_TRANS_FAIL triggers |
| `Q_BACKTRANS_MS` | Back-transition pitch limit ramp time | Controls pitch-up rate when transitioning FW→VTOL |
| `Q_TRANS_DECEL` | Deceleration for approach calculations | Used by RTL/AUTO to judge when to start VTOL transition |
| `TKOFF_THR_MIN` | Minimum throttle during transition | Applied during VTOL→FW transition |

### Assist Parameters

| Parameter | Purpose | Spearhead Notes |
|-----------|---------|-----------------|
| `Q_ASSIST_SPEED` | Airspeed below which VTOL motors assist | Set above stall speed. Use -1 to disable (suppresses pre-arm warning) |
| `Q_ASSIST_ANGLE` | Attitude error (deg) above which assist activates | Backup trigger even if above Q_ASSIST_SPEED |
| `Q_ASSIST_ALT` | AGL altitude below which assist activates | Uses rangefinder/terrain data |
| `Q_ASSIST_DELAY` | Delay before assist activates after threshold | Prevents false triggers |

**Important:** If not using an airspeed sensor, the synthetic airspeed estimate can be very inaccurate. Consider carefully whether to enable assist without a real airspeed sensor.

### RTL & Failsafe Parameters

| Parameter | Purpose |
|-----------|---------|
| `Q_RTL_MODE` | Controls hybrid RTL behavior (0–3) |
| `Q_RTL_ALT` | Altitude target for VTOL portion of hybrid RTL |
| `Q_OPTIONS` | Bitmask for many behavioral tweaks (see below) |
| `FS_SHORT_ACTN` / `FS_LONG_ACTN` | Plane failsafe actions |
| `RTL_RADIUS` | Distance at which to transition from FW to VTOL in RTL |
| `RTL_ALTITUDE` | Altitude for fixed-wing return leg |

### Weathervaning Parameters

| Parameter | Purpose |
|-----------|---------|
| `Q_WVANE_ENABLE` | Enable weathervaning (1 = nose into wind) |
| `Q_WVANE_GAIN` | Lean angle → yaw rate conversion |
| `Q_WVANE_ANG_MIN` | Minimum lean angle before weathervaning activates (default 1°) |
| `Q_FWD_THR_GAIN` | Use forward motor to help hold position in wind (v4.5+) |
| `Q_FWD_THR_USE` | Controls when forward thrust assist is used |

### Forward Motor in VTOL

| Parameter | Purpose |
|-----------|---------|
| `Q_VFWD_GAIN` | (Legacy) Use forward motor for position hold in wind |
| `Q_FWD_THR_GAIN` | (v4.5+, preferred) Forward motor position hold |
| `Q_FWD_MANTHR_MAX` | Max throttle for manual forward motor channel |
| `RCx_OPTION = 209` | Assigns an RC channel for manual forward motor in VTOL modes |

---

## 5. Q_OPTIONS Bitmask (Key Bits for Spearhead)

| Bit | Value | Effect |
|-----|-------|--------|
| 0 | +1 | Keep wings level during VTOL→FW transition (no climbing with VTOL motors) |
| 1 | +2 | Use FW takeoff for TAKEOFF command (instead of VTOL) |
| 2 | +4 | Use FW landing for LAND command |
| 3 | +8 | Respect takeoff frame altitude |
| 4 | +16 | Use fixed-wing approach for VTOL_LAND |
| 5 | +32 | Use QRTL on RC failsafe in VTOL modes |
| 7 | +128 | Force assist active at all times |
| 14 | +16384 | Use only Q_A_ANGLE_MAX for VTOL angle limits (ignore FW limits) |
| 15 | +32768 | Allow throttle descent control during VTOL auto-land |
| 16 | +65536 | Disable FW approach phase in QRTL (pure VTOL return) |
| 17 | +131072 | Allow horizontal repositioning during VTOL auto-land |
| 20 | +1048576 | Use RTL (instead of QRTL) on failsafe |

---

## 6. Flight Modes

### VTOL Modes (Quad Motor Active)

| Mode | Description | Spearhead Use |
|------|-------------|---------------|
| **QSTABILIZE** | Rate-stabilized, no altitude hold. Throttle directly controls VTOL thrust. | Basic hover testing |
| **QHOVER** | Altitude hold at mid-stick. Throttle controls climb/descent rate. | Hover testing, tuning |
| **QLOITER** | Position + altitude hold. GPS required. | Primary hover mode for ops |
| **QLAND** | Automated vertical landing | Emergency / mission end |
| **QRTL** | Return to launch in VTOL (or hybrid, depending on Q_RTL_MODE) | Failsafe / return |
| **QAUTOTUNE** | Automatic PID tuning in VTOL | Essential for initial tuning |

### Fixed Wing Modes (Recommended for QuadPlane)

| Mode | Description |
|------|-------------|
| **FBWA** | Stabilized FW flight, pilot controls attitude. **Recommended over STABILIZE** for QuadPlanes. |
| **FBWB** | Stabilized with altitude hold via pitch stick |
| **CRUISE** | GPS-stabilized cruise with altitude hold |
| **AUTO** | Autonomous mission following |
| **RTL** | Return to launch (behavior depends on Q_RTL_MODE) |
| **GUIDED** | GCS-commanded waypoints |

### Modes to AVOID in QuadPlane

**STABILIZE, ACRO, TRAINING, MANUAL** — In these modes the quad motors are disabled (except MANUAL which is intentional). The stick input is insufficient for the autopilot to determine desired attitude/climb rate for VTOL assistance. Use **FBWA** instead of STABILIZE. MANUAL immediately cuts VTOL motors — risk of stall if insufficient airspeed.

---

## 7. Transitions

### VTOL → Fixed Wing

1. Switch to any FW mode (FBWA recommended for manual transitions).
2. VTOL motors continue providing lift and stability during "transition airspeed wait" phase.
3. Forward motor (IC engine) spools up; throttle stick controls forward thrust in FBWA.
4. Elevator input during transition controls VTOL climb/descent rate (unless `Q_OPTIONS` bit 0 is set).
5. Once `AIRSPEED_MIN` is reached, VTOL motor contribution ramps down over `Q_TRANSITION_MS` (default 5 sec).
6. After ramp-down, aircraft flies as pure fixed-wing.
7. If `Q_TRANS_FAIL` is set and AIRSPEED_MIN is not reached in time, the transition aborts per `Q_TRANS_FAIL_ACT`.

**Battery warning:** During transition, all motors (VTOL + pusher) can run at very high levels simultaneously. Battery sag below 3.0V/cell (LiPo) is possible and can cause crashes. This is especially relevant for high-capacity, low-C-rating batteries. Monitor voltage during first manual transitions. Solutions: higher C-rating battery, separate VTOL and FW batteries, or use `BATT_WATT_MAX` to limit current draw.

### Fixed Wing → VTOL

1. Switch to any VTOL mode (QHOVER, QLOITER, etc.).
2. Forward motor immediately stops.
3. Control surfaces continue providing stability as the aircraft decelerates.
4. VTOL motors engage and manage attitude/altitude as speed drops.
5. When transitioning to position-hold modes (QLOITER), pitch is initially limited to 0° and relaxes to `Q_A_ANGLE_MAX` over `Q_BACKTRANS_MS` to prevent violent pitch-up and altitude gain.
6. Transitioning at high speed to QLOITER will result in nose-up pitching and altitude gain — plan accordingly.

### Typical First Flight Sequence

1. VTOL takeoff in **QLOITER** or **QHOVER**
2. Switch to **FBWA**, advance throttle above 50%, fly fixed wing
3. Switch to **QHOVER** to return to quad mode, reduce throttle to 50% for hover

---

## 8. Assisted Fixed-Wing Flight

When `Q_ASSIST_SPEED` is set to a positive value (above stall speed), the VTOL motors will automatically engage to provide lift and stability whenever:

- Airspeed drops below `Q_ASSIST_SPEED`, OR
- Attitude error exceeds `Q_ASSIST_ANGLE`, OR
- Altitude drops below `Q_ASSIST_ALT` (AGL)

Assist activates after `Q_ASSIST_DELAY` seconds. It works in all FW modes except MANUAL and ACRO.

**For initial flights**, keep `Q_ASSIST_SPEED` disabled to test basic FW functionality. Enable it later, set above stall speed.

An RC switch (`RCx_OPTION = 82`) can override assist: LOW = force off, HIGH = force on.

---

## 9. Weathervaning & Wind Hold

### Active Weathervaning

For a large-wing 25 kg aircraft, wind presents a major challenge in VTOL hover. Weathervaning automatically yaws the nose into the wind during position-controlled VTOL modes (QLOITER, QLAND, QRTL, AUTO VTOL phases — NOT QSTABILIZE or QHOVER).

Set `Q_WVANE_ENABLE = 1` (nose into wind). Start with `Q_WVANE_GAIN = 1`. If yaw oscillates, reduce gain.

**Critical:** `Q_TRIM_PITCH` must be properly set before enabling pitch-driven weathervaning (`Q_WVANE_OPTIONS` bit 0). Without correct pitch trim, you'll get continuous unwanted yawing.

QuadPlanes rarely have the same yaw authority as multirotors due to greater mass and wing area. In significant wind, expect to only be able to face into the wind.

### Forward Motor Wind Hold (v4.5+)

Set `Q_FWD_THR_GAIN` and `Q_FWD_THR_USE` to use the pusher engine to help hold position in wind during VTOL modes. This dramatically reduces the lean angle needed and the load on VTOL motors. This is the preferred method over the legacy `Q_VFWD_GAIN`.

---

## 10. Return to Launch (RTL)

### Q_RTL_MODE Options

| Mode | Behavior |
|------|----------|
| **0** | Pure FW RTL — fly back and loiter as fixed wing. VTOL motors only engage if Q_ASSIST triggers. |
| **1** | Hybrid — FW return, switch to VTOL at `RTL_RADIUS`, descend to `Q_RTL_ALT`, land vertically. |
| **2** | FW return → loiter to `Q_RTL_ALT` at `Q_FW_LND_APR_RAD` → face wind → QRTL and land. |
| **3** | FW return → approach → airbrake → VTOL transition → land. Most automated option. Also used by default for QRTL mode. |

**For Spearhead:** `Q_RTL_MODE = 3` is recommended — it provides a complete automated FW approach and VTOL landing from any distance, with airbraking to reduce speed before VTOL transition.

### QRTL Mode

By default, QRTL behaves the same as RTL with `Q_RTL_MODE = 3`. If `Q_OPTIONS` bit 16 is set, QRTL becomes pure VTOL return (no FW approach) — only use this when close to home.

### RC Failsafe in VTOL Modes

Regardless of `FS_SHORT_ACTN`/`FS_LONG_ACTN`, RC failsafe during VTOL flight triggers QLAND, QRTL, or RTL depending on `Q_OPTIONS` bits 5 and 20. If failsafe occurs during VTOL takeoff, it immediately switches to QLAND.

---

## 11. AUTO Missions

### VTOL Takeoff (NAV_VTOL_TAKEOFF)

- Set altitude parameter to desired height above takeoff point
- Aircraft climbs at `Q_WP_SPD_UP` until altitude reached
- Then proceeds to next waypoint, transitioning to FW as needed
- Latitude/longitude of the command are ignored
- Set `Q_NAVALT_MIN` to a non-zero value to prevent GPS-noise-induced roll/pitch during initial climb (forces level until that altitude)

### VTOL Landing (NAV_VTOL_LAND)

**Firmware 4.1+:** By default, NAV_VTOL_LAND remains in FW mode, flies to near the landing point, executes an airbrake maneuver, transitions to VTOL, navigates precisely to the landing point, and descends. This allows placing the land point at any distance from the last waypoint.

For a mission-end landing, place the last waypoint at ~20 m AGL before the NAV_VTOL_LAND command. For small QuadPlanes, 60–80 m separation is good; for Spearhead's size and speed, plan for a larger distance.

### Mixing FW and VTOL in Missions

Use `DO_VTOL_TRANSITION` with parameter 3 (VTOL mode) or 4 (FW mode) to switch mid-mission. Example: VTOL takeoff → FW cruise to survey area → transition to VTOL for hover photography → FW return → VTOL land.

### Hovering in Missions

Set `Q_GUIDED_MODE = 1` to make LOITER commands execute as VTOL hover instead of FW circle. This allows pausing at waypoints for photography while flying the rest of the mission as FW.

### GUIDED Mode

With `Q_GUIDED_MODE = 1`, GUIDED mode destination hold is done as VTOL hover. The approach is FW; transition to VTOL begins at `WP_LOITER_RAD` (80 m is a good starting value).

---

## 12. Setup Tips (Spearhead-Specific)

### Level Calibration

Perform accelerometer "level" calibration with the aircraft in its normal **cruise attitude** (i.e., with the wing at its cruise angle of attack, not flat). This sets the AHRS reference for fixed-wing flight.

For the VTOL motors, it's ideal to mechanically tilt each motor to be truly vertical when the airframe is in this cruise attitude. This typically requires 3–5° of forward tilt on the motor mounts. Benefits:
- Eliminates the need for `Q_TRIM_PITCH` correction
- Wings generate lift during hover in wind, reducing VTOL motor load
- Prevents sudden pitch changes during FW→VTOL transitions

If using 3D-printed motor mounts, building in this tilt is easy.

### Hover Throttle Calibration

If the aircraft doesn't hover at mid-stick (±6%) in QSTABILIZE, adjust `Q_M_THST_HOVER` to the correct percentage. Enable `Q_M_HOVER_LEARN` in QLOITER/QHOVER for automatic learning.

### Large Vehicle Considerations

- Calibrate the accelerometer on the bench before installation (Large Vehicle MagCal recommended)
- At 25 kg MTOW, the higher inertia means no benefit from raising `SCHED_LOOP_RATE` above 300 Hz
- The pitch and roll limits in VTOL modes are the lesser of `Q_A_ANGLE_MAX` or the FW limits (`PTCH_LIM_MAX_DEG`, `PTCH_LIM_MIN_DEG`, `ROLL_LIMIT_DEG`)
- `Q_BCK_PIT_LIM` limits backward pitch when airspeed is at AIRSPEED_MIN to prevent structural loads during deceleration

### QLOITER Ground Behavior

Landing detection in QLOITER is less sophisticated than Copter. If GPS shows movement while on the ground, the aircraft may try to tip over attempting to hold position. Disarm promptly after landing, or use QLAND for automated landing.

---

## 13. Simulation (SITL)

Before flying Spearhead, test all mission profiles, transitions, and failsafe behaviors in SITL:

```bash
sim_vehicle.py -j4 -v Plane -f quadplane --console --map
```

Add `-w` to reset parameters to QuadPlane defaults. Load a sample VTOL mission:

```
wp load ../Tools/autotest/Generic-Missions/KSFO-VTOL.txt
```

Use a USB transmitter adapter with the joystick module to practice QuadPlane controls.

---

## 14. Pre-Flight Checklist (Derived from ArduPilot Docs)

1. **Motor alignment verified** on bench (prop-tip distance method)
2. **Motor rotation direction confirmed** (X frame: inward toward fuselage)
3. **VTOL motor spin-up tested** in QSTABILIZE while secured
4. **Hover throttle calibrated** (`Q_M_THST_HOVER` or `Q_M_HOVER_LEARN`)
5. **Q_TRIM_PITCH set** (or motors mechanically tilted) — no backward drift in hover
6. **Q_ASSIST_SPEED set** above stall speed (or -1 for initial flights)
7. **Transition timeout** (`Q_TRANS_FAIL`) configured for safety
8. **RTL mode configured** (`Q_RTL_MODE`, `RTL_RADIUS`, `Q_RTL_ALT`)
9. **Failsafe actions verified** in SITL
10. **Airspeed sensor calibrated** (strongly recommended for transition accuracy)
11. **Battery voltage monitored** during first manual transition — check for excessive sag
12. **Weathervaning configured** if operating in wind

---

## 15. Key Documentation Links

| Topic | URL |
|-------|-----|
| QuadPlane Overview | https://ardupilot.org/plane/docs/quadplane-overview.html |
| Building | https://ardupilot.org/plane/docs/quadplane-building.html |
| Frame Setup | https://ardupilot.org/plane/docs/quadplane-frame-setup.html |
| Flying | https://ardupilot.org/plane/docs/quadplane-flying.html |
| Transitions | https://ardupilot.org/plane/docs/quadplane-transitions.html |
| Assisted Flight | https://ardupilot.org/plane/docs/assisted_fixed_wing_flight.html |
| RTL | https://ardupilot.org/plane/docs/quadplane_rtl.html |
| Parameters | https://ardupilot.org/plane/docs/quadplane-parameters.html |
| Flight Modes | https://ardupilot.org/plane/docs/quadplane-flight-modes.html |
| Weathervaning | https://ardupilot.org/plane/docs/quadplane-weathervaning.html |
| AUTO Missions | https://ardupilot.org/plane/docs/quadplane-auto-mode.html |
| Tips | https://ardupilot.org/plane/docs/quadplane-tips.html |
| Simulation | https://ardupilot.org/plane/docs/quadplane-simulation.html |
| ESC Calibration | https://ardupilot.org/plane/docs/quadplane-esc-calibration.html |
| Porter OctoQuadPlane Build | https://diydrones.com/profiles/blogs/building-flying-and-not-crashing-a-large-octaquadplane |
| Mozzie Build | http://mozzie.readthedocs.io/ |
