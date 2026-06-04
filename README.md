# GFI Toolbox — Geomorphic Flood Index 2.0

A QGIS plugin for rapid, DEM-based flood-prone area delineation using the **Geomorphic Flood Index (GFI)**.

[![QGIS](https://img.shields.io/badge/QGIS-≥3.22-green?logo=qgis)](https://plugins.qgis.org/)
[![License](https://img.shields.io/badge/License-GPL--2.0-blue)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-orange)](CHANGELOG.md)

---

## Overview

The Geomorphic Flood Index (Samela et al. 2017) maps flood-prone areas purely from terrain data — no hydraulic simulation required. This plugin implements GFI 1.0 and the hierarchical backwater extension GFI 2.0 (Manfreda et al. 2026) as a QGIS Processing tool.

MATLAB version: https://zenodo.org/records/18903835

**GFI formula:**

```
GFI = ln(hr / H)
```

Where:
- `hr = ( (A_river · Δx²) / 10⁶ )^n` — reference water depth at the mapped channel pixel
- `H` = height of the terrain pixel above its nearest channel pixel (HAND)
- A pixel is flood-prone when `GFI ≥ optimal threshold`

GFI 2.0 updates `hr` along each channel reach by propagating backwater influence from downstream confluences, improving accuracy in wide floodplains.

---

## Features

- Automatic channel network extraction with downstream propagation (Giannoni 2005)
- Strahler stream order computed internally — no GRASS or external tools required
- Hillslope-to-river and river-to-confluence topological mapping
- GFI 1.0 and GFI 2.0 computation
- ROC calibration following the `ROCcurve_maggiore` method (see MATLAB code), with configurable threshold sweep step
- Binary flood maps from each version's optimal threshold
- Optional HAND, Water Depth, Marginal Area, and channel raster outputs
- CSV performance report: AUC, CSI, Kappa, F1, TPR, TNR, FPR, FNR, Bias
- Supports GRASS r.watershed (1–8) and ESRI D8 (powers-of-2) flow direction encodings
- Save outputs to folder or as temporary QGIS layers

---

## Requirements

| Dependency | Version |
|-----------|---------|
| QGIS | ≥ 3.40 |
| Python | ≥ 3.9 (bundled with QGIS) |
| NumPy | bundled with QGIS |
| GDAL | bundled with QGIS |

No additional Python packages are required.

---

## Installation

### From the QGIS Plugin Repository (coming soon)

1. Open QGIS → **Plugins** → **Manage and Install Plugins…**
2. Search for `Geomorphic Flood Index`
3. Click **Install Plugin**

### Manual installation

1. Download the latest release ZIP from the [Releases](https://github.com/koke-research/gfi_toolbox/releases) page
2. Open QGIS → **Plugins** → **Manage and Install Plugins…** → **Install from ZIP**
3. Select the downloaded file and click **Install Plugin**

### Developer install

```bash
# Clone the repository
git clone https://github.com/koke-research/gfi_toolbox.git

# Copy (or symlink) to your QGIS plugin folder
# Windows:
xcopy gfi_toolbox "%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\gfi_toolbox" /E /I

# Linux / macOS:
cp -r gfi_toolbox ~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/
```

Then restart QGIS and enable the plugin under **Plugins** → **Manage and Install Plugins…**.

---

## Input Data

| Layer | Format | Description |
|-------|--------|-------------|
| DEM | Raster | Digital Elevation Model, any projection |
| Flow Accumulation | Raster | From GRASS `r.watershed` or equivalent |
| Flow Direction | Raster | D8 — GRASS (1–8) or ESRI (1–128) auto-detected |
| Slope | Raster | In degrees |
| Flood Reference Map | Raster (optional) | Binary (0 = dry, 1 = flooded) — required for ROC calibration |

All rasters must share the same extent, resolution, and CRS.

**Recommended preprocessing in QGIS (GRASS):**

```
r.watershed  elevation=DEM  accumulation=FlowAcc  drainage=FlowDir  ...
r.slope.aspect  elevation=DEM  slope=Slope  format=degrees
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Channel threshold (th) | 100 000 | Channel initiation: `A·S^1.7 > th` |
| Exponent (n) | 0.3544 | Scaling exponent for `hr` |
| GFI 2.0 iterations | 6 | Backwater propagation iterations |
| ROC step size | 0.01 | Threshold sweep resolution (normalised [-1,1]) |

---

## Outputs

### Always produced
| File | Description |
|------|-------------|
| `GFI_GFI_v1.tif` | GFI 1.0 index raster |
| `GFI_GFI_v2.tif` | GFI 2.0 index raster |
| `GFI_Flood_Map_v1.tif` | Binary flood map — GFI 1.0 ≥ optimal threshold |
| `GFI_Flood_Map_v2.tif` | Binary flood map — GFI 2.0 ≥ optimal threshold |

### Optional
| File | Description | Requires flood map |
|------|-------------|--------------------|
| `GFI_WD_v1.tif` | Water depth GFI 1.0 | Yes |
| `GFI_WD_v2.tif` | Water depth GFI 2.0 | Yes |
| `GFI_HAND_v1.tif` | Height Above Nearest Drainage (v1) | No |
| `GFI_HAND_v2.tif` | Height Above Nearest Drainage (v2) | No |
| `GFI_Marginal_Area.tif` | Flood reference projected onto hillslopes | Yes |
| `GFI_Channel_Network.tif` | Binary channel mask | No |
| `GFI_Strahler_Order.tif` | Strahler stream order | No |
| `GFI_Hr_v1.tif` | Reference water depth at channel | No |
| `GFI_Hr_v2.tif` | Updated reference water depth | No |
| `GFI_performance_report.csv` | Performance metrics | Yes |

### Performance report columns
`version, TP, TN, FP, FN, TPR_sensitivity, TNR_specificity, FPR, FNR, CSI_threat_score, F1_score, Kappa, Bias, AUC, a_calibrated, opt_threshold_orig, opt_threshold_norm, opt_FPR, opt_TPR, threshold_th, exponent_n, v2_iterations, roc_step_size, cellsize, run_datetime`

---

## Scientific Background

The plugin implements the methodology described in:

- **Manfreda, S., Navarro, J. S., Albertini, C., Zhuang, R., Pacia, F. D., Chaturvedi, S., & Samela, C.** (2026). Geomorphic flood index 2.0: Enhanced tools for delineating flood-prone areas in data-scarce regions. *CATENA*, 271, 110242. https://doi.org/10.1016/j.catena.2026.110242

- **Samela, C., Troy, T. J., & Manfreda, S.** (2017). Geomorphic classifiers for flood-prone areas delineation for data-scarce environments. *Advances in Water Resources*, 102, 13–28. https://doi.org/10.1016/j.advwatres.2017.01.007

- **Manfreda, S., & Samela, C.** (2019). A digital elevation model based method for a rapid estimation of flood inundation depth. *Journal of Flood Risk Management*, 12(S1), e12541. https://doi.org/10.1111/jfr3.12541

- **Giannoni F, Roth G, Rudari R.** (2005). A procedure for drainage network identification from geomorphology and its application to the prediction of the hydrologic response. *Advances in Water Resources* 28: 567 - 581. DOI: 10.1016/j.advwatres.2004.11.013

- **Navarro, J. S., Albertini, C., zhuang,  ruodan, Chaturvedi, S., Pacia, F. D., Samela, C., & Manfreda, S.** (2026). Geomorphic Flood Index (GFI) version 2.0 (Version 2.0). *Zenodo*. https://doi.org/10.5281/ZENODO.18903835


---

## Repository Structure

```
gfi_toolbox/
├── __init__.py          # QGIS plugin entry point
├── gfi_main.py          # Plugin GUI and execution pipeline
├── gfi_engine.py        # Core hydrological computations
├── metadata.txt         # QGIS plugin metadata
├── logo.png             # Plugin icon
├── LICENSE              # GPL-2.0
├── README.md
└── CHANGELOG.md
```

---

## Contributing

Bug reports and pull requests are welcome at https://github.com/jsaavedran/gfi_toolbox/issues.

Please follow the existing code style: pure NumPy, no external dependencies, and English-only strings and comments.

---

## License

This plugin is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License v2** as published by the Free Software Foundation.

See [LICENSE](LICENSE) for the full text.

---

## Authors

- Jorge Saavedra Navarro (QGIS port and GFI 2.0 implementation)
- Sadashiv Chaturvedi
- Cinzia Albertini
- Caterina Samela
- Salvatore Manfreda

Contact: nandres049@gmail.com
