# Spearhead web review

Review URL on the ranch network / ranch VPN: `http://10.3.10.123:4195/`

- Orbit, pan, zoom; 3D/top/front/side views and Fit.
- Five system visibility switches and Solo isolation.
- Click a part to see its part number, description and quantity; wireframe toggle.
- STEP and GLB downloads; mobile parts drawer.

The GLB is a tessellated preview of the same PT1 assembly, not a replacement for
the exact STEP model. CAD-to-glTF changes Z-up millimetres to Y-up metres. Mesh
bounding dimensions are checked in the browser. All geometry stays in place.

## Build and serve

From this directory:

```sh
npm ci
../.venv/bin/python export_viewer.py
../.venv/bin/python serve.py --port 4195
```

The server binds the requested port strictly on all interfaces; only this
viewer directory and the explicit generated-aircraft STEP route are exposed.
No CDN, API key or cloud upload is needed. `assets/spearhead-v3.glb` and
`node_modules/` are generated/ignored, while the exporter and package lock are
tracked. The STEP download requires `../output/spearhead.step` from the CAD
build. The live detached server logs to `/tmp/spearhead-viewer.log`.

## Verify

`node verify.mjs` uses a separate headless Chrome instance, not the operator's
browser tabs. It checks the delivered LAN URL, geometry scale, rendering
budget, subsystem controls, picking, camera motion, actual STEP download,
mobile behavior and browser errors. Evidence goes to `/tmp/spearhead-*`.

Set `VIEWER_URL`, `CHROME_PATH` and `VIEWER_EVIDENCE` to override the checked URL,
the Chrome executable (default: macOS path) and the evidence folder.
