# PT1 assembly review

Source: Fusion design `AR-SH-PT1-0000 PT1` version 35, exported without the hidden
`REF` group.

## Deliverable

- [Assembly code and instructions](../../src/design/README.md)
- [Offline preview](spearhead.png)
- [Viewer on desktop](web-viewer-desktop.png) and [on mobile](web-viewer-mobile.png)
- [Machine-readable validation](validation.json)
- [Viewer checks](web-viewer-checks.json)
- [Side review](../../src/design/pt1/pt1_side_review.csv), computed on Fusion v32; the
  tip-rib pair was re-checked on v35 and is still an exact mirror

## What the checks cover

`validate.py` compares the STEP with a body dump read from Fusion v35 and with the
part list: file hashes, 5 systems, validity of 490 solids, 426 part occurrences, 215 part
numbers, quantities per side, 482 numbered bodies, adaptive volume, centroid and bounds
of every occurrence against Fusion, the 3950 mm span, and a STEP export and import
roundtrip. `viewer/verify.mjs` checks the viewer in a headless browser.

The checks do **not** cover physical fit, interferences, joint motion, aircraft
mass or flight readiness.

## Known issues

Seven solids read as invalid in OCP from the Fusion STEP, though each is a valid solid
in Fusion. Their volumes and centroids agree with Fusion within the limits. The validity
check in `validate.py` reports them as failing on purpose; it is not loosened, so the
script exits non-zero until they are dealt with. AP214 is used because an AP242 rewrite
of the export read more invalid solids (13 against 9 on Fusion v34). OCP's automatic
shape repair fixes none of them.

| Part | Fusion volume, mm³ | Volume from STEP, mm³ | Difference | Error reported by OCP |
| --- | --- | --- | --- | --- |
| 13106L | 291830 | 291822 | -0.003 % | none named on a sub-shape |
| 13106R | 291412 | 291406 | -0.002 % | none named on a sub-shape |
| 13407L | 269933 | 269896 | -0.014 % | self-intersecting wire, unorientable face |
| 13407R | 269933 | 269896 | -0.014 % | self-intersecting wire, unorientable face |
| 14102 | 99941 | 99941 | +0.000 % | self-intersecting wire, unorientable face |
| 14103 | 140603 | 140603 | +0.000 % | self-intersecting wire, unorientable face |
| 32101L | 162259 | 162307 | +0.030 % | unorientable face |

These parts have 7 to 20 sliver faces or very short edges in Fusion, some with
zero-length edges. Writing the assembly back to STEP through Build123D changes which
solids read as invalid (12 solids after the roundtrip); `validation.json` lists both sets.

### Tip ribs, fixed in Fusion v35

On Fusion v34 the tip ribs 22208 and 23208 read 1.89 % too large from the STEP. A sliver
face of 4.0 x 0.0001 mm joined a thin wall around the leading-edge stick notch to the rest
of the rib, and the outer loop of one large face read back as two wires. Fusion v35 deletes
that sliver and removes the wall, 3.56 mm³ on each rib, at the end of the timeline. Each
rib lost 0.006 % of its volume. The stick and the leading-edge sheet did not move.

## Earlier reconstruction

This branch began as Thomas's Build123D reconstruction of the aircraft from the
separate Fusion files as they stood on September 26, 2026. It built the first complete
aircraft in one model, and the loader, checks and viewer here come from it. The
September 26 snapshot itself is no longer in the tree; it remains in this branch's
history until the pull request is squashed, and the source designs remain in Fusion.

What that work established:

- The project had no maintained complete-aircraft master, only separate designs.
- Parts needed for the aircraft were hidden in the sources: a fourth motor mount,
  six centre-wing ribs and eleven rear landing-gear bodies.
- Only one outer wing was detailed; the other had to be mirrored by inference.
- Default OCP volume integration overstated a spline leading-edge part by about
  18.8 %, against 254,117.8 mm³ with adaptive integration. The volume function here
  is the fix found then.
- A saved-visible export silently loses hidden geometry, so hidden bodies have to be
  counted.

Gaps it recorded or left, and their state in PT1:

| Reconstruction | PT1 |
| --- | --- |
| Outer wing mirrored by inference | Both wings are present, numbered and compared |
| Rear spar web hidden in the wing source, absent | 22102 and 23102 |
| Left tail spar panel absent | 31103, both sides |
| Front gear from the older `lg` v3 design | Current front gear, 15101 to 15105 |
| Propulsion and covering skins not fitted | 4 lift motors, 4 ESCs and the skins are included |
| No part numbers | 215 part numbers with materials and quantities |
