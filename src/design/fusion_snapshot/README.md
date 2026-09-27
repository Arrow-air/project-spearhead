# Fusion snapshot — 2026-09-26

All 26 current Project Spearhead design files, at their latest available cloud
revision when captured. STEP reflects saved visibility; native F3D archives
preserve the design content, including hidden geometry. Sources were not modified.

CAD assets use Git LFS (`git lfs pull`). See [assembly instructions](../README.md) for selection, inferred geometry and limitations.

| Design | Version | Saved UTC | STEP | Native |
| --- | --- | --- | --- | --- |
| conceptual_design | 14 | 2026-05-22 | [STEP](conceptual_design/source.step) | [F3D](conceptual_design/source.f3d) |
| cs_hinge | 3 | 2026-06-04 | [STEP](cs_hinge/source.step) | [F3D](cs_hinge/source.f3d) |
| DLE55RA | 8 | 2026-06-07 | [STEP](DLE55RA/source.step) | [F3D](DLE55RA/source.f3d) |
| fuselage | 24 | 2026-06-22 | [STEP](fuselage/source.step) | [F3D](fuselage/source.f3d) |
| fuselage_conceptual | 12 | 2026-06-11 | [STEP](fuselage_conceptual/source.step) | [F3D](fuselage_conceptual/source.f3d) |
| fuselage_exterior | 1 | 2026-06-12 | [STEP](fuselage_exterior/source.step) | [F3D](fuselage_exterior/source.f3d) |
| fuselage_v2 | 41 | 2026-09-22 | [STEP](fuselage_v2/source.step) | [F3D](fuselage_v2/source.f3d) |
| fuselage_v2_4mm | 1 | 2026-07-11 | [STEP](fuselage_v2_4mm/source.step) | [F3D](fuselage_v2_4mm/source.f3d) |
| KDEXF-UAS55AHVC | 1 | 2026-07-27 | [STEP](KDEXF-UAS55AHVC/source.step) | [F3D](KDEXF-UAS55AHVC/source.f3d) |
| latch | 1 | 2026-06-26 | [STEP](latch/source.step) | [F3D](latch/source.f3d) |
| le_mold_mainwing | 2 | 2026-08-04 | [STEP](le_mold_mainwing/source.step) | [F3D](le_mold_mainwing/source.f3d) |
| le_mold_tailwing | 1 | 2026-08-03 | [STEP](le_mold_tailwing/source.step) | [F3D](le_mold_tailwing/source.f3d) |
| lg | 3 | 2026-06-24 | [STEP](lg/source.step) | [F3D](lg/source.f3d) |
| mad-v8013-ipe | 1 | 2026-07-15 | [STEP](mad-v8013-ipe/source.step) | [F3D](mad-v8013-ipe/source.f3d) |
| mad_ampx_80a_5-14s | 1 | 2026-07-15 | [STEP](mad_ampx_80a_5-14s/source.step) | [F3D](mad_ampx_80a_5-14s/source.f3d) |
| MAIN SKIDS VTOL | 3 | 2026-06-05 | [STEP](MAIN_SKIDS_VTOL/source.step) | [F3D](MAIN_SKIDS_VTOL/source.f3d) |
| main_wing_structure | 42 | 2026-08-07 | [STEP](main_wing_structure/source.step) | [F3D](main_wing_structure/source.f3d) |
| MOTOR MOUNT 19x19 | 3 | 2026-06-05 | [STEP](MOTOR_MOUNT_19x19/source.step) | [F3D](MOTOR_MOUNT_19x19/source.f3d) |
| motor_mount | 7 | 2026-08-06 | [STEP](motor_mount/source.step) | [F3D](motor_mount/source.f3d) |
| stinger35cc | 2 | 2026-04-19 | [STEP](stinger35cc/source.step) | [F3D](stinger35cc/source.f3d) |
| tail_wing_structure | 2 | 2026-07-10 | [STEP](tail_wing_structure/source.step) | [F3D](tail_wing_structure/source.f3d) |
| tail_wing_structure_v2_USER_EDITS_NOW | 15 | 2026-08-05 | [STEP](tail_wing_structure_v2_USER_EDITS_NOW/source.step) | [F3D](tail_wing_structure_v2_USER_EDITS_NOW/source.f3d) |
| tattu_12s_16000mah | 2 | 2026-04-20 | [STEP](tattu_12s_16000mah/source.step) | [F3D](tattu_12s_16000mah/source.f3d) |
| temp | 1 | 2026-04-21 | [STEP](temp/source.step) | [F3D](temp/source.f3d) |
| wing_conceptual | 11 | 2026-04-22 | [STEP](wing_conceptual/source.step) | [F3D](wing_conceptual/source.f3d) |
| wing_connection | 1 | 2026-06-28 | [STEP](wing_connection/source.step) | [F3D](wing_connection/source.f3d) |

Extra recovery: [actual hidden fourth motor mount](motor_mount/hidden_fourth_mount.step).

The drawing and CSV catalog entries are not 3D design exports. Separate old,
wiring and Longshot projects are outside this snapshot.

Metadata: [catalog](catalog.json), [export index](index.json),
[native archive status](archives.json), and each design’s `manifest.json`.
The manifests distinguish native component bodies, occurrence transforms,
visibility and linked revisions. Native transform translations are centimetres;
body bounds and exported geometry are millimetres.

Reproducibility scripts in `scripts/` include Fusion-only capture scripts and
normal Python import/measurement scripts. Assembly builds do not require Fusion.
