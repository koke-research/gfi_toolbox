# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Version numbering follows [Semantic Versioning](https://semver.org/).

---

## [1.0.1] — 2026-06-13

### Fixed
- PyQt5/PyQt6 compatibility: replaced `QFrame.HLine`, `QFrame.Sunken`, `QgsFileWidget.GetDirectory`, and `QgsMapLayerProxyModel.RasterLayer` with a `try/except` block that resolves the correct enum namespace for QGIS 3.x (PyQt5) and QGIS 4.x (PyQt6)
- Replaced `QDialog.exec_()` with `exec()` — required in PyQt6 / QGIS 4.x
- `qgisMinimumVersion` corrected to `3.22`; added `qgisMaximumVersion=4.99` to allow installation on QGIS 4.x

---

## [1.0.0] — 2026-06-04

### Added
- GFI 1.0 computation: `GFI = ln(hr / H)` (Samela et al. 2017)
- GFI 2.0 computation: hierarchical backwater propagation to confluences (Albertini et al. 2022)
- Continuous channel network extraction from CI threshold `A · S^1.7 > th` with downstream gap-filling (Giannoni 2005)
- Strahler stream order computed internally — no GRASS or external tools required
- Hillslope-to-river mapping: each hillslope pixel traced downstream to its nearest channel pixel
- River-to-confluence mapping: each channel pixel traced to its next Strahler-order increase
- HAND (Height Above Nearest Drainage) rasters for both versions
- ROC calibration following the `ROCcurve_maggiore` MATLAB method:
  - GFI normalised to [-1, 1]
  - Threshold sweep from -1 to +1 with configurable step size
  - Optimisation by minimising `FPR + FNR` (balanced error)
  - AUC via rectangular averaging (`areaundercurve.m` method)
  - Calibration coefficient `a = exp(-t_original)`
- Binary flood maps from each version's optimal threshold (`GFI ≥ t`)
- Water depth outputs: `WD = max(0, hr · a − H)`
- Marginal Area raster: flood reference projected onto hillslope pixels via channel mapping
- CSV performance report with: TP, TN, FP, FN, TPR, TNR, FPR, FNR, CSI, F1, Kappa, Bias, AUC, calibration parameters
- Auto-detection of GRASS r.watershed (1–8) and ESRI D8 (powers-of-2) flow direction encodings
- Fixed `raster_to_array`: replaced raw `QByteArray` byte cast with `provider.block(band, extent, cols, rows)` + `block.value(r, c)` — eliminates garbage values (3e+348, 0)
- Configurable output mode: save to folder or temporary QGIS layers
- Scrollable output layer checklist with essential/optional distinction
- Progress bar with per-step status messages

### Plugin info
- QGIS minimum version: 3.22
- No external Python dependencies (NumPy and GDAL bundled with QGIS)
