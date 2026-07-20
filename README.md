# Project Spearhead

**An open-source fixed-wing VTOL aircraft for long-endurance missions.**

Project Spearhead is a ~25 kg MTOW hybrid quadplane developed by
[Arrow Air](https://arrowair.com): electric motors provide vertical takeoff
and landing, while an internal-combustion pusher engine is intended to provide
efficient long-range cruise.

Spearhead targets long-endurance work such as survey, rural surveillance, and
light cargo. Just as importantly, it is Arrow's first fixed-wing aircraft of
the current era — a deliberate learning platform for fixed-wing sizing, hybrid
propulsion, transition flight, and long-range operations on the path toward
larger aircraft.

## Status

**Phase 1 (electric flight validation) — PT1 prototype build in progress.**

The first prototype is being built in Ankara, Türkiye: laser-cut plywood ribs
and spars, carbon-fiber structure, 3D-printed skin and structural connectors,
and a fully connectorized electrical harness. PT1 flies electric-only; the
internal-combustion cruise system comes in Phase 2.

The aircraft is still being validated, so requirements and component choices
are provisional unless captured in a dated document in `docs/`.

## Documentation

Start with the [project intro](docs/index.md). The `docs/` folder is the
living engineering record and is published to the Arrow docs site.

| Resource | Description |
|----------|-------------|
| [Electrical Master Document](docs/electrical-master.md) | Full electrical architecture (SPH-E-001) |
| [PT1 Parts List](docs/pt1-parts-list.md) | Electrical procurement for the PT1 build |
| [Information Notes](docs/information-note/) | Aerodynamics and flight-dynamics design notes |
| [Reference Guides](docs/reference-guide/) | Stable technical references (e.g. ArduPilot QuadPlane) |

## Repository structure

```
project-spearhead/
├── docs/                        # Documentation (published to arrowair.com)
│   ├── information-note/        # Dated engineering information notes
│   └── reference-guide/         # Stable technical reference material
└── src/
    └── tools/                   # Engineering tools
        ├── initial_sizing/      # Aircraft sizing studies
        ├── propulsion_model/    # Propulsion modeling
        └── flight_dynamics_model/  # Flight dynamics / stability analysis
```

## Project phases

1. **Electric flight validation** — build and fly the airframe on electric
   VTOL propulsion only; validate hover, transition surfaces, and structure.
2. **Hybrid integration** — add the IC pusher engine, fuel and cooling
   systems; validate vibration handling, transition, and cruise flight.
3. **BVLOS and long-range testing** — telemetry range, lost-link behavior,
   and recovery protocols.
4. **Operational tailoring** — mission configurations (survey/cargo pods),
   operating limits, and a performance manual.

## Getting involved

- Read the [project intro](docs/index.md) and [CONTRIBUTING.md](CONTRIBUTING.md).
- Join the [Arrow Discord](https://discord.gg/arrow) — Spearhead has a weekly
  technical call and active build discussion.
- The original project proposal lives on the
  [Arrow DAO forum](https://dao.arrowair.com/t/project-spearhead-proposal-discussion/153).
