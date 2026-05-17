# Project Spearhead

Project Spearhead is Arrow's 25 kg MTOW hybrid fixed-wing VTOL UAV
project. It is being developed as a quadplane-style platform: electric motors
provide vertical takeoff and landing capability, while an internal-combustion
pusher engine is intended to provide efficient cruise flight.

Spearhead is aimed at long-endurance rural surveillance and limited payload or
emergency cargo missions. Just as importantly, it is a hardware learning
platform for Arrow: the project is intended to build real experience with
fixed-wing aircraft sizing, hybrid propulsion, transition flight, long-range
telemetry, safety modes, and larger aircraft systems.

These docs are a living engineering record. Spearhead is still in active
prototype development, so details will change as the team builds, tests, and
learns. Treat the proposal as the project charter and these docs as the current
state of execution.

## Project goals

Spearhead is intended to develop and validate:

- Aircraft sizing and structural integrity for a 25 kg class fixed-wing VTOL
  aircraft.
- Hybrid propulsion using electric VTOL motors for takeoff and landing, plus an
  internal-combustion engine for cruise.
- Generator integration, cruise battery-charging strategy, and energy
  management.
- Hover-to-cruise and cruise-to-landing transition control.
- Internal-combustion engine vibration mitigation and avionics protection.
- Long-range telemetry, including investigation of RF and possible
  satellite-based links.
- Modular payload support for surveillance and limited cargo use cases.
- Safety and failure handling for lost-link, engine-out, and long-range
  operations.

The goal is not to produce a polished commercial aircraft in the first pass.
The goal is to build a capable prototype platform and the internal know-how
needed for future larger systems.

## Current development phase

Spearhead is currently in early prototype development. The near-term emphasis
is electric flight validation and prototype build preparation before full
hybrid integration.

Current work includes:

- Airframe design and manufacturability.
- Wing, tail, and control-surface sizing.
- Prototype construction methods using balsa, Oratex, printed ribs, and
  laser-cut spars/ribs.
- ArduPilot QuadPlane setup, SITL, and FlightGear visualization.
- Off-the-shelf-first electrical architecture for the initial prototype.
- Servo, CAN-to-PWM, and wiring validation.
- Internal-combustion engine, starter, generator, and propulsion trade studies.
- Ground and flight-test planning around the Ankara build/test effort.

Because the aircraft is still being validated, requirements and component
choices should be read as provisional unless they are captured in a dated
decision document.

## Roadmap summary

The proposal defines four broad project phases:

### Phase 1: Electric flight validation

Build and test the full airframe with electric VTOL propulsion only. This phase
validates hover stability, control-surface authority, structural integrity, and
the wing detachment system before adding the complexity of the
internal-combustion cruise system.

### Phase 2: Hybrid integration

Install the internal-combustion pusher engine, fuel system, cooling system, and
related hybrid hardware. This phase focuses on static engine runs, temperature
and vibration measurement, flight-controller filtering, VTOL-cruise transition,
cruise flight, and battery sizing.

### Phase 3: BVLOS and long-range testing

Expand from short-range testing toward long-range flight profiles. This phase
investigates telemetry range, link latency, command reliability, lost-link
behavior, recovery protocols, and obstacle or environment awareness.

### Phase 4: Operational tailoring

Refine the platform around specific mission configurations, including cargo and
surveillance pods. This phase should produce clearer operating limits,
payload/range tradeoffs, and a performance manual while transferring Spearhead
lessons into future larger aircraft work.

## How to read these docs

These docs are organized around project execution rather than a finished
aircraft manual:

- **Status** pages describe current progress, blockers, and near-term
  priorities.
- **Roadmap** pages translate the proposal phases into the current execution
  plan.
- **Decisions** record dated technical decisions, rationale, and confidence
  level.
- **Open questions** list unresolved design, test, and operations questions.
- **Meeting notes** preserve summarized project history.
- **Research** pages collect working notes, trade studies, and experiments.
- **Reference guides** collect stable technical reference material used by the
  project.

When there is a conflict, current status and dated decision pages should be
treated as more up to date than older meeting notes or research artifacts.

## Contributing

Useful contributions right now include:

- Turning meeting notes into clear decisions, open questions, and action items.
- Reviewing aerodynamic sizing, control-surface assumptions, and transition
  risks.
- Improving ArduPilot QuadPlane, SITL, and FlightGear setup notes.
- Researching propulsion, generator, telemetry, and payload-system options.
- Helping define ground-test, taxi-test, and first-flight checklists.
- Making prototype manufacturing steps easier to repeat and document.

If you are not sure where to start, look for open questions and current status
notes first. Spearhead is moving quickly, and small clarifying contributions
can save the team a lot of rework.

## Source charter

The project started from the Arrow DAO Spearhead proposal discussion:

<https://dao.arrowair.com/t/project-spearhead-proposal-discussion/153>
