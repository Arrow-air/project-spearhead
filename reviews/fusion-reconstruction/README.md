# Fusion → Build123D reconstruction review

Date: 2026-09-26. Base commit: `e19da0485eec22891dc9027828713faa40bc6732`.
Branch: `hex/build123d-fusion-aircraft`.

## Deliverable

- [Assembly code and instructions](../../src/design/README.md)
- [Generated STEP reopened in Fusion](fusion-roundtrip.png)
- [Offline preview](spearhead.png)
- [Machine-readable validation](validation.json)
- [Cloud catalog](../../src/design/fusion_snapshot/catalog.json)
- [STEP export index](../../src/design/fusion_snapshot/index.json)
- [Native archive results](../../src/design/fusion_snapshot/archives.json)

Run `src/design/.venv/bin/python src/design/spearhead.py` from the repository
root to recreate `src/design/output/spearhead.step`.

The model contains 318 solids: 113 fuselage, 61 original wing, 55 tail,
11 boom/mount, 11 landing gear, 47 inferred mirrored wing and 3 restored hidden
motor-mount solids, plus six fixed center-wing ribs and eleven rear gear bodies.
It preserves the shared Fusion coordinates. This is a
STEP-based assembly with a deliberate mirrored-wing completion, not a claim
that the original parts have been rebuilt parametrically.

## Evidence and limitations

The current cloud project has separate structural files rather than a single
maintained aircraft master. The chosen versions and design-selection rationale
are in the assembly README. Every f3d design in that project was captured; the
separate old project and wiring project are excluded.

Native source documents remained unmodified. The fourth mount was recovered
from hidden source bodies through an unsaved disposable document. The missing
opposite outer-wing detail is explicitly inferred; verify handed features
before using it for manufacture. Current flight hardware and covering skins
are not fitted.

**Validation: 49/49 geometry checks and 22/22 browser/download checks passed.** All 318 solids survived the combined
STEP export/import roundtrip and remained valid. The generated STEP was also
reopened in native Fusion; counts are recorded in fusion-roundtrip.json.
That unsaved review document remains open; source documents were not changed. All 26 design files imported;
the separate historical fuselage has one invalid solid, excluded from the aircraft.

The checks cover asset hashes, import inventories, source-body correspondence,
solid validity, bounds, adaptive volumes, reflection selection, hidden-mount
recovery and combined STEP roundtrip. They do **not** cover physical fit,
interferences, joint motion, aircraft mass/CG or flight readiness.

### Resolved diagnostic: spline volume integration

Default non-adaptive OCP mass properties overstated a wing leading-edge body's
volume by roughly 18.8%. Native Fusion high-accuracy volume and native mesh
agreed. Adaptive OCP integration returned 254,117.8 mm³ rather than the default
301,871.2 mm³. After using adaptive integration throughout, source comparison
errors fall below 0.018% for the three principal structures.

Evidence: `fusion_snapshot/*/high-accuracy.json`, `import-metrics.json`
(includes default and adaptive values), `body-map.json`,
`native-mesh-checks.json` and `mesh-comparison.json`. Small mesh-volume
mismatches are tessellation-dependent and are not used as an exact geometry
identity test.

### Prior source issues kept separate

Historical `fuselage` v24 has an invalid imported solid; it is preserved in the
snapshot but is not selected for this assembly. A zero-volume object in the
current tail spar is omitted by STEP. `temp` was a simulation study: Design
workspace activation resolved its initial export failure.

### Review correction: center-wing ribs and rear landing gear

The initial 301-solid model omitted hidden `wing_inner` ribs and the rear
`lg` bodies from fuselage_v2 v41. Recovered six ribs and eleven rear gear
bodies via disposable DirectDesign documents without modifying the source.
Duplicate hidden spars and alternate front gear were deliberately excluded.
Exact selections, native bounds and high-accuracy volumes are preserved in
`fusion_snapshot/fuselage_v2/center-gear-recovery.json`. The corrected web
asset is `spearhead-v2.glb`; the original review URL remains unchanged.
