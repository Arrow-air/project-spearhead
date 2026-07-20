---
title: Intro
sidebar_label: Intro
sidebar_position: 1
description: Introduction to Project Spearhead
---

# Project Spearhead

Project Spearhead is Arrow's ~25 kg MTOW hybrid fixed-wing VTOL aircraft
project. It is a quadplane-style platform: electric motors provide vertical
takeoff and landing, while an internal-combustion pusher engine is intended to
provide efficient cruise flight.

Spearhead is aimed at long-endurance missions such as survey, rural
surveillance, and light cargo. Just as importantly, it is a hardware learning
platform for Arrow: the project builds real experience with fixed-wing
aircraft sizing, hybrid propulsion, transition flight, long-range telemetry,
and larger aircraft systems — a deliberate step on Arrow's path toward bigger
aircraft.

These docs are a living engineering record. Spearhead is in active prototype
development, so details will change as the team builds, tests, and learns.
Treat the [DAO proposal](https://dao.arrowair.com/t/project-spearhead-proposal-discussion/153)
as the project charter and these docs as the current state of execution.

## Project goals

Spearhead is intended to develop and validate:

- Aircraft sizing and structural integrity for a 25 kg class fixed-wing VTOL
  aircraft.
- Hybrid propulsion using electric VTOL motors plus an internal-combustion
  pusher engine for cruise.
- Generator integration, cruise battery-charging strategy, and energy
  management.
- Hover-to-cruise and cruise-to-landing transition control.
- Internal-combustion engine vibration mitigation and avionics protection.
- Long-range telemetry, including RF and possible satellite-based links.
- Modular payload support for surveillance and limited cargo use cases.
- Safety and failure handling for lost-link, engine-out, and long-range
  operations.

The goal is not a polished commercial aircraft in the first pass. The goal is
a capable prototype platform and the internal know-how needed for future
larger systems.

## Where the project is now

**Phase 1: the PT1 prototype build is underway in Ankara, Türkiye**
(as of July 2026). PT1 is the electric-only validation aircraft — the
internal-combustion cruise system comes later.

Current work includes:

- Airframe construction: laser-cut plywood ribs and spars, carbon-fiber
  structure, and 3D-printed skin and structural connectors.
- Electrical harness assembly per the
  [electrical master document](./electrical-master.md), with a fully
  connectorized, bench-tested-first approach.
- ArduPilot QuadPlane configuration on a Pixhawk flight controller, plus
  SITL simulation work.
- Aerodynamic and flight-dynamics analysis (see the
  [information notes](./information-note/)).

Because the aircraft is still being validated, requirements and component
choices are provisional unless captured in a dated document.

## Roadmap

The proposal defines four broad phases:

### Phase 1: Electric flight validation

Build and test the full airframe with electric VTOL propulsion only. This
validates hover stability, control-surface authority, structural integrity,
and the wing detachment system before adding the complexity of the
internal-combustion cruise system.

### Phase 2: Hybrid integration

Install the internal-combustion pusher engine, fuel system, cooling, and
related hybrid hardware. Focus areas: static engine runs, temperature and
vibration measurement, flight-controller filtering, VTOL–cruise transition,
cruise flight, and battery sizing.

### Phase 3: BVLOS and long-range testing

Expand from short-range testing toward long-range flight profiles: telemetry
range, link latency, command reliability, lost-link behavior, recovery
protocols, and environment awareness.

### Phase 4: Operational tailoring

Refine the platform around specific mission configurations, including cargo
and surveillance pods. This phase should produce clearer operating limits,
payload/range tradeoffs, and a performance manual while transferring
Spearhead lessons into future larger aircraft work.

## How these docs are organized

- **Engineering documents** such as the
  [electrical master document](./electrical-master.md) carry a doc ID,
  revision, author, and date. They are the authoritative specification for
  their subsystem.
- **Information notes** (`information-note/`) are dated design and analysis
  notes, organized by phase and discipline (aerodynamics, flight dynamics).
- **Reference guides** (`reference-guide/`) collect stable technical
  reference material, such as the
  [ArduPilot QuadPlane reference](./reference-guide/ardupilot-quadplane.md).
- **Planning notes** (such as the PT1 electrical planning note and
  [parts list](./pt1-parts-list.md)) capture build-cycle execution detail.

When documents conflict, the newer dated document wins.

## Contributing

Useful contributions right now include:

- Reviewing aerodynamic sizing, stability analysis, and transition risks.
- Improving ArduPilot QuadPlane, SITL, and simulation setup notes.
- Researching propulsion, generator, telemetry, and payload-system options
  for the hybrid phases.
- Helping define ground-test, taxi-test, and first-flight checklists.
- Turning meeting discussion into clear dated decisions and documentation.

See the [contributing guide](https://www.arrowair.com/docs/contributing/intro)
and join the Arrow Discord — Spearhead has a weekly technical call and active
build discussion.
