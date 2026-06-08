# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['../nondimit.py'],
    pathex=[],
    binaries=[],
    datas=[('README.md', '.'), ('aero_frame_notation.png', '.')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='Nondimit',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='Nondimit',
)
app = BUNDLE(
    coll,
    name='Nondimit.app',
    icon=None,
    bundle_identifier='com.arrowair.project-spearhead.nondimit',
)
