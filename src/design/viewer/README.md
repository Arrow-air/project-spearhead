# Spearhead web review

Review URL on the ranch network / ranch VPN: **http://10.3.10.123:4195/**

- Orbit, pan, zoom; 3D/top/front/side views and Fit.
- Seven subsystem visibility switches and Solo isolation.
- Saved-visible-only view (251 solids) versus reconstructed view (301).
- Click a part to identify it; wireframe toggle.
- STEP and GLB downloads; mobile parts drawer.
- Blue marks inferred wing detail; teal marks the recovered hidden mount.

The 18 MB GLB is a tessellated preview of the same Build123D assembly, not a
replacement for the exact STEP model. CAD-to-glTF changes Z-up millimetres to
Y-up metres. Mesh bounding dimensions are checked in the browser. Fusion's
face primitives are consolidated into 299 render meshes; the recovered mount
is one render mesh containing three solids. All geometry stays in place.

## Build and serve

From this directory:

```sh
npm ci
../.venv/bin/python export_viewer.py
../.venv/bin/python serve.py --port 4195
```

The server binds the requested port strictly on all interfaces; only this
viewer directory and the explicit generated-aircraft STEP route are exposed.
No CDN, API key or cloud upload is needed. `assets/spearhead.glb` and
`node_modules/` are generated/ignored, while the exporter and package lock are
tracked. The STEP download requires `../output/spearhead.step` from the CAD
build. The live detached server logs to `/tmp/spearhead-viewer.log`.

## Verify

`node verify.mjs` uses a separate headless Chrome instance, not the operator's
browser tabs. It checks the delivered LAN URL, geometry scale, rendering
budget, subsystem controls, picking, camera motion, actual STEP download,
mobile behavior and browser errors. Evidence goes to `/tmp/spearhead-viewer-*`.
Set the Chrome executable path in the script if running outside macOS.
