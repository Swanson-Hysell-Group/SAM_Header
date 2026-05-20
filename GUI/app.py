from __future__ import annotations

import ctypes
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from PySide6.QtCore import QObject, QThread, Qt, Signal, Slot, QTimer, QUrl
from PySide6.QtGui import QDesktopServices, QFont, QIcon
from PySide6.QtWidgets import (
    QApplication,
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QLayout,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QListWidget,
    QListWidgetItem,
    QScrollArea,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sam_generator import (
    INVALID_INPUT_DATA_MESSAGE,
    CorrectionOptions,
    ConversionResult,
    convert_csv_to_block_outputs,
    convert_csv_to_sam,
    write_stitched_block_orientation_file,
)

# --- Block sample review/correction logic ---
def correct_block_orientations(sample_df, reverse_sun=False, add_90=False, apply_declination=False, declination=0.0):
    """
    Returns a new DataFrame with corrected block strike orientations.
    """
    df = sample_df.copy()
    out = []
    for idx, row in df.iterrows():
        try:
            sun = float(row.get('shadow_angle', 'nan'))
            mag = float(row.get('magnetic_core_strike', 'nan'))
            dec = float(row.get('IGRF_local_dec', declination))
        except Exception:
            sun = mag = dec = float('nan')
        orig = sun
        if reverse_sun and not pd.isna(sun):
            sun = -sun
        if add_90 and not pd.isna(sun):
            sun = sun + 90
        if apply_declination and not pd.isna(sun) and not pd.isna(dec):
            sun = sun + dec
        # Normalize
        if not pd.isna(sun):
            sun = sun % 360
        out.append({
            'sample_name': row.get('sample_name', ''),
            'original_sun_angle': orig,
            'corrected_block_strike': sun if not pd.isna(sun) else '',
            'magnetic_core_strike': mag if not pd.isna(mag) else '',
            'declination': dec if not pd.isna(dec) else '',
        })
    return pd.DataFrame(out)

APP_ICON = ROOT / 'GUI' / 'assets' / 'sam_header_icon.svg'


def _format_pair(row: pd.Series, strike_key: str, dip_key: str) -> str:
    strike = '' if pd.isna(row.get(strike_key)) else str(row.get(strike_key)).strip()
    dip = '' if pd.isna(row.get(dip_key)) else str(row.get(dip_key)).strip()
    if strike and dip:
        return f'{strike}/{dip}'
    return strike or dip or ''


def _format_sun_summary(row: pd.Series) -> str:
    parts: list[str] = []
    shadow = '' if pd.isna(row.get('shadow_angle')) else str(row.get('shadow_angle')).strip()
    gmt = '' if pd.isna(row.get('GMT_offset')) else str(row.get('GMT_offset')).strip()
    if shadow:
        parts.append(f'shadow {shadow}')
    if gmt:
        parts.append(f'GMT {gmt}')

    date_parts = []
    for key in ['year', 'month', 'days']:
        value = '' if pd.isna(row.get(key)) else str(row.get(key)).strip()
        if value:
            date_parts.append(value)
    time_parts = []
    for key in ['hours', 'minutes']:
        value = '' if pd.isna(row.get(key)) else str(row.get(key)).strip()
        if value:
            time_parts.append(value.zfill(2) if value.isdigit() else value)

    if date_parts:
        date_text = '-'.join(date_parts)
        if time_parts:
            date_text = f"{date_text} {' : '.join(time_parts).replace(' : ', ':')}"
        parts.append(date_text)

    return ' | '.join(parts)


def _format_quality_flags(row: pd.Series) -> str:
    flags: list[str] = []
    if pd.isna(row.get('sample_name')) or str(row.get('sample_name')).strip() == '':
        flags.append('missing name')
    if pd.isna(row.get('magnetic_core_strike')) or str(row.get('magnetic_core_strike')).strip() == '':
        flags.append('no mag strike')
    if pd.isna(row.get('core_dip')) or str(row.get('core_dip')).strip() == '':
        flags.append('no core dip')
    if pd.isna(row.get('shadow_angle')) or str(row.get('shadow_angle')).strip() == '':
        flags.append('no sun')
    if pd.isna(row.get('GMT_offset')) or str(row.get('GMT_offset')).strip() == '':
        flags.append('no GMT')
    return ', '.join(flags) if flags else 'ok'


def _build_compact_orientation_preview(sample_df: pd.DataFrame) -> tuple[list[str], list[dict[str, str]]]:
    compact_columns = ['sample', 'level', 'core', 'bedding', 'sun/time', 'flags']
    compact_rows: list[dict[str, str]] = []
    for _, row in sample_df.iterrows():
        compact_rows.append(
            {
                'sample': '' if pd.isna(row.get('sample_name')) else str(row.get('sample_name')).strip(),
                'level': '' if pd.isna(row.get('strat_level')) else str(row.get('strat_level')).strip(),
                'core': _format_pair(row, 'magnetic_core_strike', 'core_dip'),
                'bedding': _format_pair(row, 'bedding_strike', 'bedding_dip'),
                'sun/time': _format_sun_summary(row),
                'flags': _format_quality_flags(row),
            }
        )
    return compact_columns, compact_rows


@dataclass(slots=True)
class CsvPreview:
    file_path: Path
    metadata: list[tuple[str, str]]
    orientation_columns: list[str]
    orientation_rows: list[dict[str, str]]
    options: CorrectionOptions
    parse_error: str | None = None


class BatchConversionWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    log_message = Signal(str)

    def __init__(self, previews: list[CsvPreview], output_directory: str) -> None:
        super().__init__()
        self.previews = previews
        self.output_directory = output_directory

    @Slot()
    def run(self) -> None:
        try:
            results: list[ConversionResult] = []
            failures: list[tuple[str, str]] = []
            stitched_rows: list[dict[str, object]] = []
            total = len(self.previews)

            for index, preview in enumerate(self.previews, start=1):
                self.log_message.emit(f'[{index}/{total}] Processing {preview.file_path.name}')
                try:
                    if preview.options.orientation_mode == 'block':
                        result = convert_csv_to_block_outputs(
                            str(preview.file_path),
                            output_directory=self.output_directory or None,
                            options=preview.options,
                            logger=self.log_message.emit,
                        )
                        stitched_rows.extend(result.stitched_rows)
                    else:
                        result = convert_csv_to_sam(
                            str(preview.file_path),
                            output_directory=self.output_directory or None,
                            options=preview.options,
                            logger=self.log_message.emit,
                        )
                except Exception:
                    failures.append((str(preview.file_path), traceback.format_exc()))
                else:
                    results.append(result)

            if stitched_rows and results:
                stitched_output = write_stitched_block_orientation_file(
                    results[0].output_directory,
                    stitched_rows,
                    logger=self.log_message.emit,
                )
                results[0].generated_files.append(stitched_output)
        except Exception:
            self.failed.emit(traceback.format_exc())
        else:
            self.finished.emit({'results': results, 'failures': failures})


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle('SAM Header Builder')
        if APP_ICON.exists():
            self.setWindowIcon(QIcon(str(APP_ICON)))
        self._set_startup_geometry()
        self._startup_geometry_applied = False

        self.worker_thread: QThread | None = None
        self.worker: BatchConversionWorker | None = None
        self.previews: list[CsvPreview] = []
        self.current_index: int = -1
        self.syncing_options = False

        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        root_layout = QVBoxLayout(central_widget)
        self.root_layout = root_layout
        root_layout.setSizeConstraint(QLayout.SizeConstraint.SetNoConstraint)
        root_layout.setContentsMargins(16, 14, 16, 14)
        root_layout.setSpacing(8)

        hero = self._build_hero()
        self.hero = hero
        root_layout.addWidget(hero)

        top_row_container = QWidget()
        self.top_row_container = top_row_container
        top_layout = QHBoxLayout(top_row_container)
        top_layout.setContentsMargins(0, 0, 0, 0)
        top_layout.setSpacing(10)
        root_layout.addWidget(top_row_container)

        inputs_panel = self._build_inputs_panel()
        preview_panel = self._build_preview_panel()
        options_panel = self._build_options_panel()
        self.inputs_panel = inputs_panel
        self.preview_panel = preview_panel
        self.options_panel = options_panel
        top_layout.addWidget(inputs_panel, 4)
        top_layout.addWidget(preview_panel, 3)
        top_layout.addWidget(options_panel, 4)

        root_layout.addSpacing(14)

        table_panel = self._build_table_panel()
        self.table_panel = table_panel
        root_layout.addWidget(table_panel, 2)

        log_panel = self._build_log_panel()
        self.log_panel = log_panel
        root_layout.addWidget(log_panel)

        footer = QLabel('Output files are always written with Windows CRLF line endings on every platform.')
        footer.setObjectName('FooterLabel')
        self.footer = footer
        root_layout.addWidget(footer)

        self._apply_styles()
        self._sync_mode_hint()

    def showEvent(self, event) -> None:  # type: ignore[override]
        super().showEvent(event)
        if not self._startup_geometry_applied:
            self._startup_geometry_applied = True
            # Apply one more pass after layout + frame metrics are known.
            QTimer.singleShot(0, self._fit_window_within_screen)

    def _set_startup_geometry(self) -> None:
        screen = QApplication.primaryScreen()
        if screen is None:
            self.resize(1040, 700)
            return

        available = screen.availableGeometry()

        width = min(1248, available.width() - 24)
        height = min(760, available.height() - 36)

        width = max(700, min(width, available.width()))
        height = max(500, min(height, available.height()))

        x = available.x() + max(0, (available.width() - width) // 2)
        y = available.y() + max(0, (available.height() - height) // 2)
        self.setGeometry(x, y, width, height)

    def _fit_window_within_screen(self) -> None:
        screen = self.screen() or QApplication.primaryScreen()
        if screen is None:
            return

        available = screen.availableGeometry()
        frame = self.frameGeometry()

        target_frame_w = min(max(frame.width(), int(available.width() * 0.78)), max(200, available.width() - 12))
        target_frame_h = max(200, int(available.height() * 0.78))

        frame_margin_w = max(0, frame.width() - self.width())
        frame_margin_h = max(0, frame.height() - self.height())

        needed_client_h = self._startup_content_height_hint()
        target_client_h = max(320, max(target_frame_h - frame_margin_h, needed_client_h))

        target_client_w = max(700, target_frame_w - frame_margin_w)

        target_client_w = min(target_client_w, available.width())
        target_client_h = min(target_client_h, max(320, available.height() - frame_margin_h - 12))

        # Temporarily cap size so layout minimums cannot immediately grow the window.
        self.setMinimumSize(0, 0)
        self.setMaximumSize(target_client_w, target_client_h)
        self.setFixedHeight(target_client_h)
        self.resize(target_client_w, target_client_h)

        frame_after = self.frameGeometry()
        x = available.x() + max(0, (available.width() - frame_after.width()) // 2)
        y = available.y() + max(0, (available.height() - frame_after.height()) // 2)
        self.move(x, y)

        # Restore user-resizable behavior after startup sizing is applied.
        QTimer.singleShot(120, self._release_startup_size_cap)

    def _release_startup_size_cap(self) -> None:
        self.setMinimumHeight(400)
        self.setMaximumSize(16777215, 16777215)

    def _startup_content_height_hint(self) -> int:
        margins = self.root_layout.contentsMargins()
        spacing = self.root_layout.spacing()
        return (
            margins.top()
            + margins.bottom()
            + self.hero.sizeHint().height()
            + self.top_row_container.sizeHint().height()
            + self.table_panel.minimumHeight()
            + self.log_panel.minimumHeight()
            + self.footer.sizeHint().height()
            + (spacing * 4)
            + 12
        )

    def _wrap_panel_scroll(self, panel: QWidget) -> QScrollArea:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        panel.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        scroll.setWidget(panel)
        return scroll

    def _build_hero(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName('HeroCard')
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.setSpacing(4)

        title = QLabel('SAM Header Builder')
        title.setObjectName('HeroTitle')
        title_font = QFont('Avenir Next', 26)
        title_font.setBold(True)
        title.setFont(title_font)

        subtitle = QLabel(
            'Wrap the existing SAM generator in a desktop workflow with explicit declination controls '
            'for core samples, block samples, and bedding.'
        )
        subtitle.setWordWrap(True)
        subtitle.setObjectName('HeroSubtitle')

        layout.addWidget(title)
        layout.addWidget(subtitle)
        return frame

    def _build_inputs_panel(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName('PanelCard')
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(8)

        header = QLabel('Inputs')
        header.setObjectName('SectionTitle')
        layout.addWidget(header)

        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(8)
        grid.setColumnStretch(1, 1)   # line-edit column absorbs all flexible space
        grid.setColumnMinimumWidth(2, 130)  # button column never shrinks below text

        self.input_edit = QLineEdit()
        self.input_edit.setReadOnly(True)
        self.output_edit = QLineEdit()

        self.mode_combo = QComboBox()
        self.mode_combo.addItem('Core sample workflow', 'core')
        self.mode_combo.addItem('Block sample corrected CSV workflow', 'block')
        self.mode_combo.currentIndexChanged.connect(self._update_current_options)
        self.mode_combo.currentIndexChanged.connect(self._sync_mode_hint)

        browse_input = QPushButton('Browse CSV Files')
        browse_input.clicked.connect(self._browse_input)
        browse_output = QPushButton('Browse Folder')
        browse_output.clicked.connect(self._browse_output)

        grid.addWidget(QLabel('CSV files'), 0, 0)
        grid.addWidget(self.input_edit, 0, 1)
        grid.addWidget(browse_input, 0, 2)
        grid.addWidget(QLabel('Output folder'), 1, 0)
        grid.addWidget(self.output_edit, 1, 1)
        grid.addWidget(browse_output, 1, 2)
        grid.addWidget(QLabel('Workflow'), 2, 0)
        grid.addWidget(self.mode_combo, 2, 1, 1, 2)
        layout.addLayout(grid)

        toggle_header = QLabel('Loaded files')
        toggle_header.setObjectName('SectionTitle')
        layout.addWidget(toggle_header)

        file_hint = QLabel('Single-click to preview. Ctrl/Cmd-click for bulk settings.')
        file_hint.setObjectName('HintLabel')
        file_hint.setWordWrap(True)
        layout.addWidget(file_hint)

        self.file_list = QListWidget()
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.file_list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.file_list.setMinimumHeight(48)
        self.file_list.setMaximumHeight(72)
        self.file_list.currentRowChanged.connect(self._switch_preview)
        layout.addWidget(self.file_list)

        button_row = QHBoxLayout()
        button_row.setSpacing(6)
        self.generate_button = QPushButton('Generate Current File')
        self.generate_button.clicked.connect(self._start_current_conversion)
        self.generate_batch_button = QPushButton('Generate Batch')
        self.generate_batch_button.clicked.connect(self._start_batch_conversion)
        self.open_folder_button = QPushButton('Open Output Folder')
        self.open_folder_button.clicked.connect(self._open_output_folder)
        self.open_folder_button.setEnabled(False)
        button_row.addWidget(self.generate_button)
        button_row.addWidget(self.generate_batch_button)
        button_row.addWidget(self.open_folder_button)
        button_row.addStretch(1)
        layout.addLayout(button_row)
        return frame

    def _build_options_panel(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName('PanelCard')
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(8)

        header = QLabel('Declination Options')
        header.setObjectName('SectionTitle')
        layout.addWidget(header)

        self.prefer_sun_checkbox = QCheckBox('Prefer sun compass data when available')
        self.prefer_sun_checkbox.setChecked(True)
        self.prefer_sun_checkbox.toggled.connect(self._update_current_options)
        self.correct_core_checkbox = QCheckBox('Apply local declination when magnetic orientation becomes core strike')
        self.correct_core_checkbox.setChecked(True)
        self.correct_core_checkbox.toggled.connect(self._update_current_options)
        self.correct_block_checkbox = QCheckBox('Apply local declination when magnetic orientation becomes block strike')
        self.correct_block_checkbox.toggled.connect(self._update_current_options)
        self.correct_bedding_checkbox = QCheckBox('Allow bedding strike correction from local declination')
        self.correct_bedding_checkbox.setChecked(True)
        self.correct_bedding_checkbox.toggled.connect(self._update_current_options)
        self.reverse_sun_checkbox = QCheckBox('Reverse sun compass reading (block)')
        self.reverse_sun_checkbox.setChecked(False)
        self.reverse_sun_checkbox.toggled.connect(self._update_current_options)
        self.add_90_checkbox = QCheckBox('Add 90° to orientation')
        self.add_90_checkbox.setChecked(False)
        self.add_90_checkbox.toggled.connect(self._update_current_options)

        layout.addWidget(self.prefer_sun_checkbox)
        layout.addWidget(self.correct_core_checkbox)
        layout.addWidget(self.correct_block_checkbox)
        layout.addWidget(self.correct_bedding_checkbox)

        self.mode_hint = QLabel()
        self.mode_hint.setWordWrap(True)
        self.mode_hint.setObjectName('HintLabel')
        layout.addWidget(self.mode_hint)
        layout.addWidget(self.reverse_sun_checkbox)
        layout.addWidget(self.add_90_checkbox)

        apply_row = QHBoxLayout()
        apply_row.setSpacing(6)
        self.apply_selected_button = QPushButton('Apply To Selected')
        self.apply_selected_button.clicked.connect(self._apply_current_options_to_selected)
        self.apply_all_button = QPushButton('Apply To All')
        self.apply_all_button.clicked.connect(self._apply_current_options_to_all)
        apply_row.addWidget(self.apply_selected_button)
        apply_row.addWidget(self.apply_all_button)
        layout.addLayout(apply_row)
        layout.addStretch(1)
        return frame

    def _build_preview_panel(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName('PanelCard')
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(6)

        self.preview_title = QLabel('CSV preview')
        self.preview_title.setObjectName('SectionTitle')
        layout.addWidget(self.preview_title)

        self.preview_error = QLabel('Select CSV files to preview site metadata and validation details.')
        self.preview_error.setObjectName('HintLabel')
        self.preview_error.setWordWrap(True)
        layout.addWidget(self.preview_error)

        metadata_label = QLabel('Site metadata')
        metadata_label.setObjectName('SubSectionTitle')
        layout.addWidget(metadata_label)

        self.metadata_table = QTableWidget(0, 2)
        self.metadata_table.setHorizontalHeaderLabels(['Field', 'Value'])
        self.metadata_table.verticalHeader().setVisible(False)
        self.metadata_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.metadata_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.metadata_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.metadata_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.metadata_table.setMinimumHeight(126)
        self.metadata_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.metadata_table)
        layout.addStretch(1)
        return frame

    def _build_table_panel(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName('PanelCard')
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(6)

        orient_label = QLabel('Input Orientation Data')
        orient_label.setObjectName('SectionTitle')
        layout.addWidget(orient_label)

        table_hint = QLabel('Showing the full CSV sample table for the currently selected file.')
        table_hint.setObjectName('HintLabel')
        table_hint.setWordWrap(True)
        layout.addWidget(table_hint)

        self.orientation_table = QTableWidget(0, 0)
        self.orientation_table.verticalHeader().setVisible(False)
        self.orientation_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.orientation_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.orientation_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.orientation_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.orientation_table, 1)

        self.block_table = QTableWidget(0, 0)
        self.block_table.verticalHeader().setVisible(False)
        self.block_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.block_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.block_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.block_table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.block_table, 1)
        self.block_table.hide()
        frame.setMinimumHeight(190)
        return frame

    def _build_log_panel(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName('PanelCard')
        frame.setMinimumHeight(118)
        frame.setMaximumHeight(150)
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(4)

        header = QLabel('Run log')
        header.setObjectName('SectionTitle')
        layout.addWidget(header)

        self.log_output = QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.log_output.setMinimumHeight(72)
        layout.addWidget(self.log_output, 1)
        return frame

    def _apply_styles(self) -> None:
        self.setStyleSheet(
            '''
            QWidget {
                background: #ffffff;
                color: #293222;
                font-family: "Segoe UI", "Helvetica Neue", sans-serif;
                font-size: 13px;
            }
            QMainWindow {
                background: qlineargradient(
                    x1: 0, y1: 0, x2: 1, y2: 1,
                    stop: 0 #ffffff,
                    stop: 0.5 #eef6df,
                    stop: 1 #d9ebb8
                );
            }
            QFrame#HeroCard {
                background: qlineargradient(
                    x1: 0, y1: 0, x2: 1, y2: 1,
                    stop: 0 #4f8c3f,
                    stop: 0.55 #68a94f,
                    stop: 1 #85bf63
                );
                border-radius: 22px;
            }
            QFrame#PanelCard {
                background: rgba(255, 255, 255, 0.97);
                border: 1px solid rgba(118, 157, 80, 0.18);
                border-radius: 18px;
            }
            QLabel {
                background: transparent;
            }
            QLabel#HeroTitle {
                color: #f8fff4;
                font-size: 29px;
                letter-spacing: 1px;
            }
            QLabel#HeroSubtitle {
                color: rgba(247, 255, 241, 0.96);
                font-size: 13px;
            }
            QLabel#SectionTitle {
                color: #4d7036;
                font-size: 14px;
                font-weight: 700;
            }
            QLabel#SubSectionTitle {
                color: #5f8246;
                font-size: 12px;
                font-weight: 700;
            }
            QLabel#HintLabel, QLabel#FooterLabel {
                color: #677356;
            }
            QLineEdit, QComboBox, QListWidget {
                min-width: 0px;
            }
            QLineEdit, QComboBox, QPlainTextEdit, QTableWidget {
                background: #ffffff;
                border: 1px solid #d5dfb8;
                border-radius: 10px;
                padding: 4px 6px;
                selection-background-color: #74ad57;
                selection-color: #f7fff4;
            }
            QListWidget {
                background: #ffffff;
                border: 1px solid #d5dfb8;
                border-radius: 10px;
                padding: 3px;
                outline: none;
            }
            QListWidget::item {
                padding: 5px 7px;
                border-radius: 8px;
                color: #4d7036;
            }
            QListWidget::item:selected {
                background: #e8f0d0;
                color: #355224;
            }
            QListWidget::item:hover {
                background: #f1f7e0;
            }
            QComboBox::drop-down {
                border: none;
            }
            QPushButton {
                background: #6aa14b;
                color: #f8fff4;
                border: none;
                border-radius: 12px;
                padding: 6px 10px;
                font-weight: 600;
            }
            QPushButton:hover {
                background: #5a9140;
            }
            QPushButton:disabled {
                background: #c1ccb0;
                color: #eef4ea;
            }
            QPushButton#FileToggle {
                border-radius: 16px;
                padding: 6px 12px;
                background: #e8f0d0;
                color: #4d7036;
            }
            QPushButton#FileToggle:checked {
                background: #74ad57;
                color: #f8fff4;
            }
            QCheckBox {
                spacing: 6px;
                color: #31402a;
            }
            QCheckBox::indicator {
                width: 14px;
                height: 14px;
            }
            QCheckBox::indicator:unchecked {
                border: 1px solid #96ad7a;
                border-radius: 5px;
                background: #ffffff;
            }
            QCheckBox::indicator:checked {
                border: 1px solid #5a9140;
                border-radius: 5px;
                background: #88c45f;
            }
            QHeaderView::section {
                background: #edf4d7;
                color: #4a6f34;
                border: none;
                border-right: 1px solid #d7e3bf;
                border-bottom: 1px solid #d7e3bf;
                padding: 5px;
                font-weight: 700;
            }
            QScrollBar:vertical {
                background: transparent;
                width: 12px;
                margin: 4px 2px 4px 2px;
            }
            QScrollBar::handle:vertical {
                background: rgba(106, 161, 75, 0.42);
                min-height: 28px;
                border-radius: 6px;
            }
            QScrollBar::handle:vertical:hover {
                background: rgba(90, 145, 64, 0.62);
            }
            QScrollBar:horizontal {
                background: transparent;
                height: 12px;
                margin: 2px 4px 2px 4px;
            }
            QScrollBar::handle:horizontal {
                background: rgba(106, 161, 75, 0.42);
                min-width: 28px;
                border-radius: 6px;
            }
            QScrollBar::handle:horizontal:hover {
                background: rgba(90, 145, 64, 0.62);
            }
            QScrollBar::add-line, QScrollBar::sub-line,
            QScrollBar::add-page, QScrollBar::sub-page {
                background: transparent;
                border: none;
            }
            '''
        )

    def _append_log(self, message: str) -> None:
        self.log_output.appendPlainText(message)

    def _browse_input(self) -> None:
        file_names, _ = QFileDialog.getOpenFileNames(self, 'Choose CSV templates', str(ROOT), 'CSV files (*.csv)')
        if not file_names:
            return

        self.previews = [self._parse_csv_preview(Path(file_name).resolve()) for file_name in file_names]
        self.input_edit.setText(f'{len(self.previews)} file(s) selected')
        self._rebuild_file_toggles()

        if not self.output_edit.text().strip():
            self.output_edit.setText(str(Path(file_names[0]).resolve().parent))

        parse_failures = [preview for preview in self.previews if preview.parse_error]
        if parse_failures:
            self._append_log('Some files failed preview parsing and will be skipped in batch:')
            for preview in parse_failures:
                self._append_log(f'- {preview.file_path.name}: {preview.parse_error}')

            if any(preview.parse_error == INVALID_INPUT_DATA_MESSAGE for preview in parse_failures):
                QMessageBox.warning(self, 'Invalid CSV file', INVALID_INPUT_DATA_MESSAGE)

        if self.previews:
            self.file_list.setCurrentRow(0)

    def _browse_output(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, 'Choose output folder', self.output_edit.text().strip() or str(ROOT))
        if directory:
            self.output_edit.setText(directory)

    def _parse_csv_preview(self, file_path: Path) -> CsvPreview:
        default_options = CorrectionOptions(
            orientation_mode='core',
            apply_declination_to_core_strike=True,
            apply_declination_to_block_strike=False,
            apply_declination_to_bedding_strike=True,
            prefer_sun_compass=True,
        )
        try:
            if file_path.stat().st_size == 0:
                raise ValueError(INVALID_INPUT_DATA_MESSAGE)

            hdf = pd.read_csv(str(file_path), header=0, index_col=0, nrows=5, usecols=[0, 1], dtype=object)
            if hdf.empty:
                raise ValueError(INVALID_INPUT_DATA_MESSAGE)

            metadata = []
            for key, value in hdf.iloc[:, 0].items():
                metadata.append((str(key), '' if pd.isna(value) else str(value)))

            sample_df = pd.read_csv(str(file_path), header=6, dtype=object)
            if sample_df.empty:
                raise ValueError(INVALID_INPUT_DATA_MESSAGE)

            orientation_columns = list(sample_df.columns)
            orientation_rows = sample_df.fillna('').astype(str).to_dict(orient='records')

            return CsvPreview(
                file_path=file_path,
                metadata=metadata,
                orientation_columns=orientation_columns,
                orientation_rows=orientation_rows,
                options=default_options,
                parse_error=None,
            )
        except Exception as exc:
            return CsvPreview(
                file_path=file_path,
                metadata=[],
                orientation_columns=[],
                orientation_rows=[],
                options=default_options,
                parse_error=INVALID_INPUT_DATA_MESSAGE if isinstance(exc, ValueError) else str(exc),
            )

    def _rebuild_file_toggles(self) -> None:
        self.file_list.clear()
        for preview in self.previews:
            item = QListWidgetItem(preview.file_path.name)
            item.setToolTip(str(preview.file_path))
            self.file_list.addItem(item)

    @Slot(int)
    def _switch_preview(self, index: int) -> None:
        if index < 0 or index >= len(self.previews):
            return
        self.current_index = index

        preview = self.previews[index]
        self.preview_title.setText(f'CSV preview: {preview.file_path.name}')

        if preview.parse_error:
            self.preview_error.setText(f'Preview unavailable for this file: {preview.parse_error}')
            self.metadata_table.setRowCount(0)
            self.orientation_table.setRowCount(0)
            self.block_table.setRowCount(0)
            self.block_table.hide()
        else:
            self.preview_error.setText(
                f'Loaded {len(preview.orientation_rows)} row(s) from {preview.file_path.name}. '
                'Site metadata is summarized here and the full sample table is shown below.'
            )
            self._populate_metadata_table(preview)
            self._populate_orientation_table(preview)
            # If block mode, show block correction table
            if self.mode_combo.currentData() == 'block':
                self._populate_block_table(preview)
                self.block_table.show()
            else:
                self.block_table.hide()

        self._apply_options_to_controls(preview.options)
        self._sync_mode_hint()

    def _populate_block_table(self, preview: CsvPreview) -> None:
        # Show block correction table
        try:
            df = pd.read_csv(str(preview.file_path), header=6, dtype=object)
        except Exception:
            self.block_table.setRowCount(0)
            return
        reverse_sun = self.reverse_sun_checkbox.isChecked()
        add_90 = self.add_90_checkbox.isChecked()
        apply_decl = self.correct_block_checkbox.isChecked()
        # Use first non-blank declination value if present
        try:
            decl = float(df['IGRF_local_dec'].dropna().iloc[0])
        except Exception:
            decl = 0.0
        result = correct_block_orientations(df, reverse_sun, add_90, apply_decl, decl)
        self.block_table.setColumnCount(len(result.columns))
        self.block_table.setHorizontalHeaderLabels(list(result.columns))
        self.block_table.setRowCount(len(result))
        for row_index, row in enumerate(result.itertuples(index=False)):
            for col_index, value in enumerate(row):
                self.block_table.setItem(row_index, col_index, QTableWidgetItem(str(value)))
        self.block_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.block_table.horizontalHeader().setStretchLastSection(True)

    def _populate_metadata_table(self, preview: CsvPreview) -> None:
        self.metadata_table.setRowCount(len(preview.metadata))
        for row_index, (key, value) in enumerate(preview.metadata):
            self.metadata_table.setItem(row_index, 0, QTableWidgetItem(key))
            self.metadata_table.setItem(row_index, 1, QTableWidgetItem(value))
        self._resize_metadata_table()

    def _populate_orientation_table(self, preview: CsvPreview) -> None:
        self.orientation_table.clear()
        self.orientation_table.setColumnCount(len(preview.orientation_columns))
        self.orientation_table.setHorizontalHeaderLabels(preview.orientation_columns)
        self.orientation_table.setRowCount(len(preview.orientation_rows))
        for row_index, row in enumerate(preview.orientation_rows):
            for col_index, column_name in enumerate(preview.orientation_columns):
                item = QTableWidgetItem(row.get(column_name, ''))
                item.setToolTip(row.get(column_name, ''))
                self.orientation_table.setItem(row_index, col_index, item)

        self.orientation_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.orientation_table.horizontalHeader().setStretchLastSection(False)
        self.orientation_table.verticalHeader().setDefaultSectionSize(24)

    def _resize_metadata_table(self) -> None:
        header_height = self.metadata_table.horizontalHeader().height()
        row_height = self.metadata_table.verticalHeader().defaultSectionSize()
        frame = self.metadata_table.frameWidth() * 2
        visible_rows = max(5, self.metadata_table.rowCount())
        target_height = header_height + (row_height * visible_rows) + frame + 6
        self.metadata_table.setMinimumHeight(target_height)
        self.metadata_table.setMaximumHeight(target_height)

    def _sync_mode_hint(self) -> None:
        mode = self.mode_combo.currentData()
        if mode == 'block':
            self.mode_hint.setText(
                'Block workflow writes corrected CSV files for each loaded block file and a stitched '
                'orientation CSV across the batch. No SAM, sample, or .inp files are generated.'
            )
        else:
            self.mode_hint.setText(
                'Core workflow matches the standard SAM path. Sun compass values still win when present unless you disable that above.'
            )

    def _build_options_from_controls(self) -> CorrectionOptions:
        return CorrectionOptions(
            orientation_mode=str(self.mode_combo.currentData()),
            apply_declination_to_core_strike=bool(self.correct_core_checkbox.isChecked()),
            apply_declination_to_block_strike=bool(self.correct_block_checkbox.isChecked()),
            apply_declination_to_bedding_strike=bool(self.correct_bedding_checkbox.isChecked()),
            prefer_sun_compass=bool(self.prefer_sun_checkbox.isChecked()),
            reverse_block_sun_compass=bool(self.reverse_sun_checkbox.isChecked()),
            add_ninety_to_block_strike=bool(self.add_90_checkbox.isChecked()),
        )

    def _apply_options_to_controls(self, options: CorrectionOptions) -> None:
        self.syncing_options = True
        self.mode_combo.setCurrentIndex(0 if options.orientation_mode == 'core' else 1)
        self.prefer_sun_checkbox.setChecked(options.prefer_sun_compass)
        self.correct_core_checkbox.setChecked(options.apply_declination_to_core_strike)
        self.correct_block_checkbox.setChecked(options.apply_declination_to_block_strike)
        self.correct_bedding_checkbox.setChecked(options.apply_declination_to_bedding_strike)
        self.reverse_sun_checkbox.setChecked(options.reverse_block_sun_compass)
        self.add_90_checkbox.setChecked(options.add_ninety_to_block_strike)
        self.syncing_options = False

    def _update_current_options(self) -> None:
        if self.syncing_options:
            return
        if self.current_index < 0 or self.current_index >= len(self.previews):
            return
        self.previews[self.current_index].options = self._build_options_from_controls()

    def _selected_preview_indices(self) -> list[int]:
        indices = sorted(index.row() for index in self.file_list.selectionModel().selectedRows())
        if not indices and self.current_index >= 0:
            return [self.current_index]
        return indices

    def _apply_current_options_to_selected(self) -> None:
        self._apply_current_options_to_indices(self._selected_preview_indices(), 'selected')

    def _apply_current_options_to_all(self) -> None:
        self._apply_current_options_to_indices(list(range(len(self.previews))), 'all loaded')

    def _apply_current_options_to_indices(self, indices: list[int], scope_label: str) -> None:
        self._update_current_options()
        if not indices:
            QMessageBox.information(self, 'No files selected', 'Select one or more files in the file list first.')
            return

        options = self._build_options_from_controls()
        for index in indices:
            if 0 <= index < len(self.previews):
                self.previews[index].options = CorrectionOptions(
                    orientation_mode=options.orientation_mode,
                    apply_declination_to_core_strike=options.apply_declination_to_core_strike,
                    apply_declination_to_block_strike=options.apply_declination_to_block_strike,
                    apply_declination_to_bedding_strike=options.apply_declination_to_bedding_strike,
                    prefer_sun_compass=options.prefer_sun_compass,
                    reverse_block_sun_compass=options.reverse_block_sun_compass,
                    add_ninety_to_block_strike=options.add_ninety_to_block_strike,
                )

        self._append_log(f'Applied current declination settings to {len(indices)} {scope_label} file(s).')

    def _start_current_conversion(self) -> None:
        if self.current_index < 0 or self.current_index >= len(self.previews):
            QMessageBox.warning(self, 'Missing CSV', 'Choose at least one CSV template first.')
            return
        self._start_worker([self.previews[self.current_index]])

    def _start_batch_conversion(self) -> None:
        if not self.previews:
            QMessageBox.warning(self, 'Missing CSV', 'Choose at least one CSV template first.')
            return
        self._start_worker(self.previews)

    def _start_worker(self, previews: list[CsvPreview]) -> None:
        self._update_current_options()
        output_path = self.output_edit.text().strip()

        valid_previews = [preview for preview in previews if preview.parse_error is None]
        if not valid_previews:
            QMessageBox.warning(self, 'No valid CSV files', 'No selected CSV files could be parsed.')
            return

        self.log_output.clear()
        self._append_log(f'Starting conversion for {len(valid_previews)} file(s)...')
        self.generate_button.setEnabled(False)
        self.generate_batch_button.setEnabled(False)
        self.open_folder_button.setEnabled(False)

        self.worker_thread = QThread(self)
        self.worker = BatchConversionWorker(valid_previews, output_path)
        self.worker.moveToThread(self.worker_thread)

        self.worker_thread.started.connect(self.worker.run)
        self.worker.log_message.connect(self._append_log)
        self.worker.finished.connect(self._handle_success)
        self.worker.failed.connect(self._handle_failure)
        self.worker.finished.connect(self.worker_thread.quit)
        self.worker.failed.connect(self.worker_thread.quit)
        self.worker_thread.finished.connect(self._cleanup_worker)
        self.worker_thread.start()

    @Slot(object)
    def _handle_success(self, payload: object) -> None:
        if not isinstance(payload, dict):
            self._handle_failure('Unexpected worker payload.')
            return

        results: list[ConversionResult] = payload.get('results', [])
        failures: list[tuple[str, str]] = payload.get('failures', [])

        self._append_log('')
        generated_files = sum(len(result.generated_files) for result in results)
        self._append_log(f'Successful files: {len(results)} | Failed files: {len(failures)}')
        self._append_log(f'Generated artifacts: {generated_files}')

        if results:
            self.output_edit.setText(str(results[0].output_directory))

        if failures:
            self._append_log('')
            self._append_log('Failures:')
            for file_name, trace in failures:
                self._append_log(f'- {Path(file_name).name}')
                self._append_log(trace)

        self.generate_button.setEnabled(True)
        self.generate_batch_button.setEnabled(True)
        self.open_folder_button.setEnabled(True)
        QMessageBox.information(
            self,
            'Conversion complete',
            f'Completed {len(results)} file(s) with {len(failures)} failure(s).',
        )

    @Slot(str)
    def _handle_failure(self, error_text: str) -> None:
        self._append_log('')
        self._append_log(error_text)
        self.generate_button.setEnabled(True)
        self.generate_batch_button.setEnabled(True)
        QMessageBox.critical(self, 'Conversion failed', error_text)

    @Slot()
    def _cleanup_worker(self) -> None:
        if self.worker is not None:
            self.worker.deleteLater()
            self.worker = None
        if self.worker_thread is not None:
            self.worker_thread.deleteLater()
            self.worker_thread = None

    def _open_output_folder(self) -> None:
        directory = self.output_edit.text().strip()
        if not directory:
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(directory))


def main() -> int:
    if sys.platform.startswith('win'):
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('samheader.builder.desktop')

    app = QApplication(sys.argv)
    app.setApplicationName('SAM Header Builder')
    if APP_ICON.exists():
        app.setWindowIcon(QIcon(str(APP_ICON)))
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == '__main__':
    raise SystemExit(main())