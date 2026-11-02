# Contributing to Project Spearhead

Start with the general
[Arrow Contributing Guide](https://www.arrowair.com/docs/contributing/intro)
on the Arrow website — it covers how Arrow works as a DAO, how grants and
bounties are structured, and general contribution standards.

## Where the work happens

- **Discord** — join the [Arrow Discord](https://discord.gg/arrow). Spearhead
  has a weekly technical call plus ongoing build discussion; that is where
  design decisions are made and where you should introduce yourself.
- **GitHub** — this repository holds the documentation and engineering tools.
  Issues and pull requests are welcome.
- **DAO forum** — the original
  [project proposal](https://dao.arrowair.com/t/project-spearhead-proposal-discussion/153)
  is the project charter.

## What's useful right now

Spearhead is in Phase 1 (electric flight validation) with the PT1 prototype
build underway. High-value contributions include:

- Reviewing aerodynamic sizing, stability analysis, and transition-flight
  assumptions in the information notes.
- Improving the ArduPilot QuadPlane reference and SITL/simulation setup notes.
- Research on propulsion, generator, telemetry, and payload options for the
  hybrid phases.
- Helping define ground-test, taxi-test, and first-flight checklists.
- Turning meeting discussion into clear dated decisions and documentation.
- Improving the engineering tools in `src/tools/`.

If you're unsure where to start, ask in Discord — small, concrete
contributions with clear evidence (analysis, sources, test data) are always
preferred over broad proposals.

## Documentation conventions

- Engineering documents carry a doc ID (e.g. `SPH-E-001` for the electrical
  master document), a revision, an author, and a date.
- Requirements and component choices are provisional unless captured in a
  dated document — when you find a conflict, the newer dated document wins.
- Information notes live in `docs/information-note/`, organized by phase and
  discipline.

## Pull requests

- Keep PRs focused and small enough to review.
- Commit messages follow [Conventional Commits](https://www.conventionalcommits.org/)
  (enforced by commitlint), e.g. `docs: add PT1 harness notes`.
- Spelling is checked with cspell; add legitimate project terms to
  `.cspell.project-words.txt`.
