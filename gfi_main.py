# gfi_main.py
import os
import csv
import tempfile
import numpy as np
from datetime import datetime
from osgeo import gdal

from qgis.PyQt.QtWidgets import (
    QAction, QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QMessageBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QProgressBar, QCheckBox, QScrollArea, QWidget,
    QComboBox, QFrame,
)
from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import Qgis, QgsMapLayerProxyModel, QgsRasterLayer, QgsProject
from qgis.gui import QgsMapLayerComboBox, QgsFileWidget

# PyQt5 / PyQt6 enum compatibility
try:
    # PyQt6 (QGIS 4.x)
    _HLINE       = QFrame.Shape.HLine
    _SUNKEN      = QFrame.Shadow.Sunken
    _GET_DIR     = QgsFileWidget.StorageMode.GetDirectory
    _PROXY_RASTER = QgsMapLayerProxyModel.Filter.RasterLayer
except AttributeError:
    # PyQt5 (QGIS 3.x)
    _HLINE       = QFrame.HLine          # type: ignore[attr-defined]
    _SUNKEN      = QFrame.Sunken         # type: ignore[attr-defined]
    _GET_DIR     = _GET_DIR   # type: ignore[attr-defined]
    _PROXY_RASTER = QgsMapLayerProxyModel.RasterLayer  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Output layer registry
#   id, display label, essential (always on + locked)
#
# Essential (always produced):
#   gfi_v1, gfi_v2          — GFI index rasters
#   flood_v1, flood_v2      — Binary flood maps from optimal threshold
#
# Optional (require flood map where noted):
#   wd_v1, wd_v2            — Water depth  [needs flood map]
#   h_v1, h_v2              — HAND / height above channel
#   marg_area               — Marginal Area  [needs flood map]
#   channel                 — Channel network mask
#   strahler, hr_v1, hr_v2  — Intermediate hydrological rasters
#   row_ch, col_ch, ...     — Topology index rasters
# ---------------------------------------------------------------------------
OUTPUT_LAYERS = [
    # id              label                                           essential
    ("gfi_v1",     "GFI 1.0  [ln(hr/H)]",                           True),
    ("gfi_v2",     "GFI 2.0  [hierarchical backwater]",              True),
    ("flood_v1",   "Flood Map 1.0  [optimal threshold]",             True),
    ("flood_v2",   "Flood Map 2.0  [optimal threshold]",             True),
    ("wd_v1",      "Water Depth — WD 1.0",                           False),
    ("wd_v2",      "Water Depth — WD 2.0",                           False),
    ("h_v1",       "HAND 1.0  [height above nearest channel]",       False),
    ("h_v2",       "HAND 2.0  [height above updated channel]",       False),
    ("marg_area",  "Marginal Area  [projected flood reference]",     False),
    ("channel",    "Channel Network",                                False),
    ("strahler",   "Strahler Stream Order",                          False),
    ("hr_v1",      "River Water Depth — Hr 1.0",                     False),
    ("hr_v2",      "River Water Depth — Hr 2.0",                     False),
    ("row_ch",     "Row index — nearest channel pixel",              False),
    ("col_ch",     "Col index — nearest channel pixel",              False),
    ("row_conf",   "Row index — next confluence",                    False),
    ("col_conf",   "Col index — next confluence",                    False),
]

# Outputs that require a flood map to be provided
FLOOD_MAP_REQUIRED = {"flood_v1", "flood_v2", "wd_v1", "wd_v2", "marg_area"}


# ---------------------------------------------------------------------------
# GDAL helper
# ---------------------------------------------------------------------------
def _save_raster(array, out_path, geotransform, wkt_crs,
                 nodata=-9999.0, gdal_dtype=gdal.GDT_Float32):
    rows, cols = array.shape
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(out_path, cols, rows, 1, gdal_dtype,
                       options=['COMPRESS=LZW', 'TILED=YES'])
    ds.SetGeoTransform(geotransform)
    ds.SetProjection(wkt_crs)
    band = ds.GetRasterBand(1)
    clean = np.where(np.isnan(array.astype(np.float64)), nodata, array)
    band.WriteArray(clean.astype(np.float32))
    band.SetNoDataValue(nodata)
    ds.FlushCache()
    ds = None


# ---------------------------------------------------------------------------
# Plugin entry point
# ---------------------------------------------------------------------------
class GFIPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = os.path.dirname(__file__)
        self.action = None

    def initGui(self):
        self.action = QAction("Run GFI 2.0", self.iface.mainWindow())
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("&GFI 2.0", self.action)

    def unload(self):
        self.iface.removePluginMenu("&GFI 2.0", self.action)

    # -----------------------------------------------------------------------
    # Dialog
    # -----------------------------------------------------------------------
    def run(self):
        self.dialog = QDialog(self.iface.mainWindow())
        self.dialog.setWindowTitle("GFI 2.0 — Configuration")
        self.dialog.setMinimumWidth(580)

        root = QVBoxLayout()
        root.setSpacing(8)

        # ── Section 1: Inputs ────────────────────────────────────────────────
        grp_in = QGroupBox("1. Input Layers")
        lay_in = QGridLayout()
        lay_in.setColumnStretch(1, 1)

        def _combo():
            cb = QgsMapLayerComboBox()
            cb.setFilters(_PROXY_RASTER)
            return cb

        lay_in.addWidget(QLabel("DEM:"),                        0, 0)
        self.combo_dem   = _combo(); lay_in.addWidget(self.combo_dem,   0, 1)
        lay_in.addWidget(QLabel("Flow Accumulation:"),          1, 0)
        self.combo_acc   = _combo(); lay_in.addWidget(self.combo_acc,   1, 1)
        lay_in.addWidget(QLabel("Flow Direction (D8):"),        2, 0)
        self.combo_dir   = _combo(); lay_in.addWidget(self.combo_dir,   2, 1)
        lay_in.addWidget(QLabel("Slope (degrees):"),            3, 0)
        self.combo_slope = _combo(); lay_in.addWidget(self.combo_slope, 3, 1)

        # Flood map — optional
        sep = QFrame(); sep.setFrameShape(_HLINE); sep.setFrameShadow(_SUNKEN)
        lay_in.addWidget(sep, 4, 0, 1, 2)

        flood_label = QLabel("Flood Reference Map (optional):")
        flood_label.setToolTip(
            "Binary raster (0=dry, 1=flooded). Required for ROC calibration, "
            "WD outputs, MargArea and the performance CSV.")
        lay_in.addWidget(flood_label, 5, 0)

        flood_row = QHBoxLayout()
        self.check_flood = QCheckBox("Use flood map")
        self.check_flood.setChecked(False)
        self.check_flood.toggled.connect(self._on_flood_toggled)
        flood_row.addWidget(self.check_flood)
        self.combo_flood = _combo()
        self.combo_flood.setEnabled(False)
        flood_row.addWidget(self.combo_flood)
        flood_row.addStretch()
        flood_widget = QWidget()
        flood_widget.setLayout(flood_row)
        lay_in.addWidget(flood_widget, 5, 1)

        grp_in.setLayout(lay_in)
        root.addWidget(grp_in)

        # ── Section 2: Parameters ────────────────────────────────────────────
        grp_par = QGroupBox("2. Parameters")
        lay_par = QGridLayout()

        lay_par.addWidget(QLabel("Channel threshold (th):"), 0, 0)
        self.spin_th = QDoubleSpinBox()
        self.spin_th.setRange(100, 1e8); self.spin_th.setDecimals(0)
        self.spin_th.setValue(100_000); self.spin_th.setSuffix("  m²·S^1.7")
        lay_par.addWidget(self.spin_th, 0, 1)

        lay_par.addWidget(QLabel("Exponent (n):"), 1, 0)
        self.spin_n = QDoubleSpinBox()
        self.spin_n.setRange(0.01, 1.0); self.spin_n.setDecimals(4)
        self.spin_n.setValue(0.3544)
        lay_par.addWidget(self.spin_n, 1, 1)

        lay_par.addWidget(QLabel("GFI 2.0 iterations:"), 2, 0)
        self.spin_iter = QSpinBox()
        self.spin_iter.setRange(1, 20); self.spin_iter.setValue(6)
        lay_par.addWidget(self.spin_iter, 2, 1)

        lay_par.addWidget(QLabel("ROC step size:"), 3, 0)
        self.spin_roc_steps = QDoubleSpinBox()
        self.spin_roc_steps.setRange(0.001, 0.1)
        self.spin_roc_steps.setDecimals(3)
        self.spin_roc_steps.setValue(0.01)
        self.spin_roc_steps.setSingleStep(0.005)
        self.spin_roc_steps.setToolTip(
            "Threshold sweep step for ROC calibration (normalised [-1, 1]).\n"
            "Smaller = finer search, slower. Default 0.01.")
        lay_par.addWidget(self.spin_roc_steps, 3, 1)

        grp_par.setLayout(lay_par)
        root.addWidget(grp_par)

        # ── Section 3: Outputs ───────────────────────────────────────────────
        grp_out = QGroupBox("3. Output Layers")
        lay_out = QVBoxLayout()

        # Storage mode
        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Storage:"))
        self.combo_storage = QComboBox()
        self.combo_storage.addItem("Save to folder", "folder")
        self.combo_storage.addItem("Temporary files", "temp")
        self.combo_storage.currentIndexChanged.connect(self._on_storage_changed)
        mode_row.addWidget(self.combo_storage)
        mode_row.addStretch()
        lay_out.addLayout(mode_row)

        # Folder picker
        self.folder_widget = QWidget()
        folder_row = QHBoxLayout(); folder_row.setContentsMargins(0, 0, 0, 0)
        folder_row.addWidget(QLabel("Folder:"))
        self.folder_out = QgsFileWidget()
        self.folder_out.setStorageMode(_GET_DIR)
        folder_row.addWidget(self.folder_out)
        self.folder_widget.setLayout(folder_row)
        lay_out.addWidget(self.folder_widget)

        # CSV performance report checkbox
        self.check_csv = QCheckBox("Export performance report (CSV)")
        self.check_csv.setChecked(True)
        self.check_csv.setToolTip("Requires flood reference map. Saves AUC, CSI, Kappa, TPR, TNR, etc.")
        lay_out.addWidget(self.check_csv)

        sep2 = QFrame(); sep2.setFrameShape(_HLINE); sep2.setFrameShadow(_SUNKEN)
        lay_out.addWidget(sep2)

        # Scrollable layer checklist
        scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setMaximumHeight(240)
        scroll_inner = QWidget(); scroll_layout = QVBoxLayout(); scroll_layout.setSpacing(2)

        self.output_checks = {}
        for layer_id, label, essential in OUTPUT_LAYERS:
            cb = QCheckBox(label)
            cb.setChecked(True)
            if essential:
                cb.setEnabled(False)
                cb.setToolTip("Required output.")
            if layer_id in FLOOD_MAP_REQUIRED:
                cb.setChecked(False)
                cb.setToolTip("Requires flood reference map.")
            self.output_checks[layer_id] = cb
            scroll_layout.addWidget(cb)

        btn_row = QHBoxLayout()
        btn_all  = QPushButton("Select all");      btn_all.clicked.connect(lambda: self._set_all_checks(True))
        btn_none = QPushButton("Deselect optional"); btn_none.clicked.connect(lambda: self._set_all_checks(False))
        btn_row.addWidget(btn_all); btn_row.addWidget(btn_none); btn_row.addStretch()
        scroll_layout.addLayout(btn_row)
        scroll_inner.setLayout(scroll_layout); scroll.setWidget(scroll_inner)
        lay_out.addWidget(scroll)

        grp_out.setLayout(lay_out)
        root.addWidget(grp_out)

        # ── Progress + Run ───────────────────────────────────────────────────
        self.progress_bar = QProgressBar(); self.progress_bar.setValue(0)
        root.addWidget(self.progress_bar)

        self.btn_run = QPushButton("▶  Run GFI Modelling")
        self.btn_run.setStyleSheet(
            "background-color:#2b78e4;color:white;font-weight:bold;padding:9px;font-size:13px;")
        self.btn_run.clicked.connect(self.execute_engine)
        root.addWidget(self.btn_run)

        self.dialog.setLayout(root)
        self._on_storage_changed()
        self.dialog.exec()

    # -----------------------------------------------------------------------
    # UI helpers
    # -----------------------------------------------------------------------
    def _on_flood_toggled(self, checked):
        self.combo_flood.setEnabled(checked)
        # Enable/disable flood-map-dependent outputs
        for lid in FLOOD_MAP_REQUIRED:
            self.output_checks[lid].setEnabled(checked)
            if not checked:
                self.output_checks[lid].setChecked(False)
        self.check_csv.setEnabled(checked)
        if not checked:
            self.check_csv.setChecked(False)

    def _on_storage_changed(self):
        self.folder_widget.setVisible(self.combo_storage.currentData() == "folder")
        self.dialog.adjustSize()

    def _set_all_checks(self, state):
        for layer_id, _, essential in OUTPUT_LAYERS:
            cb = self.output_checks[layer_id]
            if not essential and cb.isEnabled():
                cb.setChecked(state)

    def update_progress(self, val, msg):
        self.progress_bar.setValue(val)
        self.iface.messageBar().pushMessage("GFI 2.0", msg, level=Qgis.Info)
        QCoreApplication.processEvents()

    # -----------------------------------------------------------------------
    # Output path helper
    # -----------------------------------------------------------------------
    def _out_path(self, suffix, out_folder):
        """Return a file path, either in out_folder or as a named temp file."""
        if out_folder:
            return os.path.join(out_folder, f"GFI_{suffix}.tif")
        tmp = tempfile.NamedTemporaryFile(suffix=f"_GFI_{suffix}.tif", delete=False)
        path = tmp.name; tmp.close()
        self._temp_files.append(path)
        return path

    # -----------------------------------------------------------------------
    # Main execution
    # -----------------------------------------------------------------------
    def execute_engine(self):
        self.btn_run.setEnabled(False)
        self._temp_files = []

        try:
            # ── 0. Collect UI inputs ─────────────────────────────────────────
            layer_dem   = self.combo_dem.currentLayer()
            layer_acc   = self.combo_acc.currentLayer()
            layer_dir   = self.combo_dir.currentLayer()
            layer_slope = self.combo_slope.currentLayer()

            if not all([layer_dem, layer_acc, layer_dir, layer_slope]):
                raise ValueError("Please select all 4 input raster layers.")

            has_flood   = self.check_flood.isChecked()
            layer_flood = self.combo_flood.currentLayer() if has_flood else None
            if has_flood and layer_flood is None:
                raise ValueError("Flood reference map is checked but no layer is selected.")

            storage_mode = self.combo_storage.currentData()
            if storage_mode == "folder":
                out_folder = self.folder_out.filePath()
                if not out_folder:
                    raise ValueError("Please select an output folder.")
                os.makedirs(out_folder, exist_ok=True)
            else:
                out_folder = None

            th_val    = self.spin_th.value()
            n_val     = self.spin_n.value()
            iter_val  = self.spin_iter.value()
            roc_steps = self.spin_roc_steps.value()
            wanted   = {lid for lid, _, _ in OUTPUT_LAYERS
                        if self.output_checks[lid].isChecked()}
            save_csv = self.check_csv.isChecked() and has_flood

            from .gfi_engine import (
                raster_to_array, validate_gfi_inputs,
                build_continuous_channel, compute_strahler_order,
                hillslope_to_river, river_to_confluence,
                compute_gfi_v1, calibrate_gfi,
                compute_gfi_v2, compute_performance,
            )

            # ── 1. Read rasters ──────────────────────────────────────────────
            self.update_progress(5,  "Reading DEM…")
            dem_array, cellsize, _ = raster_to_array(layer_dem)
            self.update_progress(10, "Reading Flow Accumulation…")
            acc_array, _, _        = raster_to_array(layer_acc)
            self.update_progress(14, "Reading Flow Direction…")
            dir_array, _, _        = raster_to_array(layer_dir)
            self.update_progress(18, "Reading Slope…")
            slope_array, _, _      = raster_to_array(layer_slope)

            flood_array = None
            if has_flood:
                self.update_progress(21, "Reading Flood Reference Map…")
                flood_raw, _, _ = raster_to_array(layer_flood)
                flood_array = (flood_raw > 0).astype(np.float32)
                flood_array[np.isnan(flood_raw)] = np.nan

            # ── 1.5. Validate ────────────────────────────────────────────────
            self.update_progress(23, "Validating inputs…")
            ok, msg = validate_gfi_inputs(
                dem_array, acc_array, dir_array, flood_array, slope_array)
            if not ok:
                raise ValueError(msg)

            provider     = layer_dem.dataProvider()
            extent       = provider.extent()
            geotransform = (extent.xMinimum(), cellsize, 0,
                            extent.yMaximum(), 0, -cellsize)
            wkt_crs      = layer_dem.crs().toWkt()

            # ── 2. Channel network ───────────────────────────────────────────
            self.update_progress(25, "Extracting channel network…")
            channel = build_continuous_channel(
                acc_array, dir_array, slope_array, cellsize, th_val)

            # ── 3. Strahler order ────────────────────────────────────────────
            self.update_progress(35, "Computing Strahler order…")
            strahler = compute_strahler_order(channel, dir_array)

            # ── 4. Hillslope-to-river ────────────────────────────────────────
            self.update_progress(44, "Hillslope → nearest channel…")
            row_ch, col_ch = hillslope_to_river(
                dem_array, channel, dir_array, cellsize)

            # ── 5. River-to-confluence ───────────────────────────────────────
            self.update_progress(53, "Channel → next confluence…")
            row_conf, col_conf = river_to_confluence(
                channel, strahler, dir_array, cellsize)

            # ── 6. GFI v1.0 ─────────────────────────────────────────────────
            self.update_progress(60, "Computing GFI v1.0…")
            gfi_v1, hr_v1, h_v1, _ = compute_gfi_v1(
                dem_array, acc_array, channel, row_ch, col_ch, cellsize, n_val)

            # ── 7. Calibration v1 (requires flood map) ───────────────────────
            a_v1      = 1.0
            thr_v1    = 0.0
            marg_area = None
            roc_v1    = None
            wd_v1     = np.full_like(gfi_v1, np.nan)

            if has_flood:
                self.update_progress(65, "Calibrating GFI 1.0 (ROC)…")
                a_v1, thr_v1, marg_area, roc_v1 = calibrate_gfi(
                    gfi_v1, flood_array, row_ch, col_ch, channel,
                    version_label="1.0", roc_steps=roc_steps)
                wd_v1 = np.maximum(0.0, hr_v1 * a_v1 - h_v1)
                wd_v1 = np.where(np.isnan(h_v1), np.nan, wd_v1).astype(np.float32)
            else:
                self.iface.messageBar().pushMessage(
                    "GFI 2.0",
                    "No flood map — flood maps and WD outputs will not be produced. "
                    "GFI v2 will use a=1.0 (uncalibrated).",
                    level=Qgis.Warning)

            # Binary flood maps from optimal threshold:
            # flood pixel where GFI >= optimal_threshold
            flood_v1 = np.where(
                np.isfinite(gfi_v1),
                (gfi_v1 >= thr_v1).astype(np.float32),
                np.nan
            ).astype(np.float32)

            # ── 8. GFI v2.0 ─────────────────────────────────────────────────
            self.update_progress(70, "Computing GFI 2.0 (backwater propagation)…")
            gfi_v2, hr_v2, h_v2, _ = compute_gfi_v2(
                dem_array, acc_array, channel, strahler,
                row_ch, col_ch, row_conf, col_conf,
                cellsize, n_val, a_v1, max_iter=iter_val)

            wd_v2  = np.full_like(gfi_v2, np.nan)
            roc_v2 = None
            a_v2   = a_v1
            thr_v2 = thr_v1

            if has_flood:
                self.update_progress(78, "Calibrating GFI 2.0 (ROC)…")
                a_v2, thr_v2, _, roc_v2 = calibrate_gfi(
                    gfi_v2, flood_array, row_ch, col_ch, channel,
                    version_label="2.0", roc_steps=roc_steps)
                wd_v2 = np.maximum(0.0, hr_v2 * a_v2 - h_v2)
                wd_v2 = np.where(np.isnan(h_v2), np.nan, wd_v2).astype(np.float32)

            flood_v2 = np.where(
                np.isfinite(gfi_v2),
                (gfi_v2 >= thr_v2).astype(np.float32),
                np.nan
            ).astype(np.float32)

            # ── 9. Performance metrics & CSV ─────────────────────────────────
            perf_rows = []
            if has_flood:
                self.update_progress(83, "Computing performance metrics…")
                # Performance evaluated on binary flood maps (GFI >= optimal threshold)
                metrics_v1 = compute_performance(
                    flood_v1, flood_array, marg_area, "1.0", roc_v1)
                metrics_v2 = compute_performance(
                    flood_v2, flood_array, marg_area, "2.0", roc_v2)

                run_info = {
                    'run_datetime':  datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    'threshold_th':  th_val,
                    'exponent_n':    n_val,
                    'v2_iterations': iter_val,
                    'roc_step_size': roc_steps,
                    'cellsize':      round(cellsize, 4),
                }
                for m in [metrics_v1, metrics_v2]:
                    is_v1 = m['version'] == '1.0'
                    row = {**run_info,
                           'a_calibrated':        a_v1   if is_v1 else a_v2,
                           'opt_threshold_orig':  thr_v1 if is_v1 else thr_v2,
                           'opt_threshold_norm':  round(roc_v1['opt_t_norm'], 4) if is_v1
                                                  else round(roc_v2['opt_t_norm'], 4),
                           'opt_FPR':             round(roc_v1['opt_fpr'], 4) if is_v1
                                                  else round(roc_v2['opt_fpr'], 4),
                           'opt_TPR':             round(roc_v1['opt_tpr'], 4) if is_v1
                                                  else round(roc_v2['opt_tpr'], 4),
                           **m}
                    perf_rows.append(row)

                if save_csv and out_folder:
                    csv_path = os.path.join(out_folder, "GFI_performance_report.csv")
                    with open(csv_path, 'w', newline='') as f:
                        writer = csv.DictWriter(f, fieldnames=list(perf_rows[0].keys()))
                        writer.writeheader()
                        writer.writerows(perf_rows)
                    self.iface.messageBar().pushMessage(
                        "GFI 2.0", f"Performance report saved: {csv_path}", level=Qgis.Info)

            # ── 10. Save & load selected rasters ─────────────────────────────
            self.update_progress(88, "Saving output rasters…")

            all_arrays = {
                "gfi_v1"    : gfi_v1,
                "gfi_v2"    : gfi_v2,
                "flood_v1"  : flood_v1,
                "flood_v2"  : flood_v2,
                "wd_v1"     : wd_v1,
                "wd_v2"     : wd_v2,
                "h_v1"      : h_v1,
                "h_v2"      : h_v2,
                "marg_area" : marg_area if marg_area is not None else np.full_like(gfi_v1, np.nan),
                "channel"   : channel.astype(np.float32),
                "strahler"  : strahler.astype(np.float32),
                "hr_v1"     : hr_v1,
                "hr_v2"     : hr_v2,
                "row_ch"    : row_ch,
                "col_ch"    : col_ch,
                "row_conf"  : row_conf,
                "col_conf"  : col_conf,
            }

            suffixes = {
                "gfi_v1"    : "GFI_v1",
                "gfi_v2"    : "GFI_v2",
                "flood_v1"  : "Flood_Map_v1",
                "flood_v2"  : "Flood_Map_v2",
                "wd_v1"     : "WD_v1",
                "wd_v2"     : "WD_v2",
                "h_v1"      : "HAND_v1",
                "h_v2"      : "HAND_v2",
                "marg_area" : "Marginal_Area",
                "channel"   : "Channel_Network",
                "strahler"  : "Strahler_Order",
                "hr_v1"     : "Hr_v1",
                "hr_v2"     : "Hr_v2",
                "row_ch"    : "Row_Channel",
                "col_ch"    : "Col_Channel",
                "row_conf"  : "Row_Confluence",
                "col_conf"  : "Col_Confluence",
            }

            display_names = {
                "gfi_v1"    : "GFI 1.0",
                "gfi_v2"    : "GFI 2.0",
                "flood_v1"  : "GFI — Flood Map 1.0",
                "flood_v2"  : "GFI — Flood Map 2.0",
                "wd_v1"     : "GFI — Water Depth 1.0",
                "wd_v2"     : "GFI — Water Depth 2.0",
                "h_v1"      : "GFI — HAND 1.0",
                "h_v2"      : "GFI — HAND 2.0",
                "marg_area" : "GFI — Marginal Area",
                "channel"   : "GFI — Channel Network",
                "strahler"  : "GFI — Strahler Order",
                "hr_v1"     : "GFI — Hr 1.0",
                "hr_v2"     : "GFI — Hr 2.0",
                "row_ch"    : "GFI — Row Channel",
                "col_ch"    : "GFI — Col Channel",
                "row_conf"  : "GFI — Row Confluence",
                "col_conf"  : "GFI — Col Confluence",
            }

            for layer_id in wanted:
                if layer_id in FLOOD_MAP_REQUIRED and not has_flood:
                    continue
                arr  = all_arrays[layer_id]
                path = self._out_path(suffixes[layer_id], out_folder)
                _save_raster(arr, path, geotransform, wkt_crs)
                rl = QgsRasterLayer(path, display_names[layer_id])
                QgsProject.instance().addMapLayer(rl)

            self.update_progress(100, "Done!")
            self.iface.messageBar().pushMessage(
                "GFI 2.0",
                f"Modelling complete — {len(wanted)} layer(s) loaded."
                + (" | Performance CSV saved." if save_csv and out_folder else ""),
                level=Qgis.Success)

        except Exception as e:
            import traceback
            QMessageBox.critical(
                self.dialog, "GFI Modelling Error",
                f"{e}\n\n{traceback.format_exc()}")
        finally:
            self.btn_run.setEnabled(True)
