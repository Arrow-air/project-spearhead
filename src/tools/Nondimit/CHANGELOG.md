# Changelog

## v1.1

- Started the next versioned working folder after the v1.0 baseline.
- Added visible application version text in the window title, header, and About dialog.
- Added correction controls in Convert Raw Table and Modify Export.
- Corrections can be applied in body or wind frame.
- Each force and moment component has editable multiplier and delta values.
- Wind-frame corrections project through alpha and beta, apply the selected values, and reconstruct body-axis force and moment vectors so body and wind outputs stay consistent.
- Correction settings are stored in the export metadata options row.
- Added macOS source launcher and macOS app build script for producing a macOS ZIP package on a Mac.
- Split platform package files into `windows/` and `macos/`, and removed the Windows batch launcher.
- Moved macOS build instructions out of the main README and into `macos/BUILD.md`.
