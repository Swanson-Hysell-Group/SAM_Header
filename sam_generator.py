from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, field
from datetime import datetime as dt
from functools import reduce
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from mk_sam_utilities import igrf, sundec, to_year_fraction


WINDOWS_NEWLINE = '\r\n'
Logger = Callable[[str], None]
INVALID_INPUT_DATA_MESSAGE = 'There is no valid input data in the loaded CSV file.'


@dataclass(slots=True)
class CorrectionOptions:
    orientation_mode: str = 'core'
    apply_declination_to_core_strike: bool = True
    apply_declination_to_block_strike: bool = False
    apply_declination_to_bedding_strike: bool = True
    prefer_sun_compass: bool = True
    reverse_block_sun_compass: bool = False
    add_ninety_to_block_strike: bool = False

    def magnetic_declination_enabled(self) -> bool:
        if self.orientation_mode == 'block':
            return self.apply_declination_to_block_strike
        return self.apply_declination_to_core_strike

    def orientation_label(self) -> str:
        if self.orientation_mode == 'block':
            return 'block strike'
        return 'core strike'


@dataclass(slots=True)
class ConversionResult:
    output_directory: Path
    site_id: str
    generated_files: list[Path] = field(default_factory=list)
    stitched_rows: list[dict[str, object]] = field(default_factory=list)


def _default_logger(message: str) -> None:
    print(message)


def _is_missing(value: object) -> bool:
    return bool(pd.isna(value))


def _to_float(value: object) -> float:
    if _is_missing(value):
        return float('nan')
    return float(value)


def _yes(value: object, *, default: bool = False) -> bool:
    if _is_missing(value):
        return default
    return str(value).strip().lower() == 'yes'


def _normalize_degrees(value: float) -> float:
    normalized = value % 360
    if normalized < 0:
        normalized += 360
    return normalized


def _stringify_csv_value(value: object) -> str:
    if _is_missing(value):
        return ''
    if isinstance(value, list):
        return str(list(value)).replace(',', ';')
    return str(value)


def _read_text_with_fallback(file_name: str) -> str:
    try:
        with open(file_name, 'r', encoding='utf-8') as csv_file:
            return csv_file.read()
    except UnicodeDecodeError:
        with open(file_name, 'r', encoding='ISO-8859-1') as csv_file:
            return csv_file.read()


def _write_text(path: Path, content: str) -> None:
    with open(path, 'w', encoding='utf-8', newline='') as target_file:
        target_file.write(content)


def fix_line_breaks(file_name: str) -> None:
    csv_str = _read_text_with_fallback(file_name)
    fixed_lines = csv_str.replace('\r\n', '\n').replace('\r', '\n')
    with open(file_name, 'w', encoding='utf-8', newline='') as new_csv_file:
        new_csv_file.write(fixed_lines)


def _load_input_frames(file_name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df_cols = [
        'sample_name', 'comment', 'strat_level',
        'magnetic_core_strike', 'core_dip', 'bedding_strike',
        'bedding_dip', 'correct_bedding_using_local_dec',
        'mass', 'runs', 'sun_core_strike', 'calculated_IGRF',
        'IGRF_local_dec', 'calculated_mag_dec', 'core_strike',
        'corrected_bedding_strike',
    ]
    sdf_cols = [
        'sample_name', 'shadow_angle', 'GMT_offset',
        'year', 'month', 'days', 'hours', 'minutes',
    ]

    hdf = pd.read_csv(file_name, header=0, index_col=0, nrows=5, usecols=[0, 1])
    df = pd.read_csv(
        file_name,
        header=6,
        index_col=0,
        usecols=df_cols,
        dtype=object,
    ).transpose()
    sdf = pd.read_csv(
        file_name,
        header=6,
        index_col=0,
        usecols=sdf_cols,
        dtype=object,
    ).transpose()

    if hdf.empty or df.empty or sdf.empty:
        raise ValueError(INVALID_INPUT_DATA_MESSAGE)

    site_info = hdf.iloc[:, 0] if not hdf.empty else pd.Series(dtype=object)
    site_id = site_info.get('site_id') if not site_info.empty else None
    if pd.isna(site_id) or str(site_id).strip() == '':
        raise ValueError(INVALID_INPUT_DATA_MESSAGE)

    return hdf, df, sdf


def _calculate_declinations(
    hdf: pd.DataFrame,
    df: pd.DataFrame,
    sdf: pd.DataFrame,
    options: CorrectionOptions,
    logger: Logger,
) -> None:
    samples = df.keys()
    time_types = ['year', 'month', 'days', 'hours', 'minutes']

    logger('---------------------LOCAL MAGNETIC DECLINATION-----------------------')

    for sample in samples:
        if not sdf[sample].isnull().any():
            time_values = [str(int(sdf[sample][time_type])) for time_type in time_types]
            assert len(time_values[0]) == 4, (
                'must input full year for sun compass calculation (i.e. YYYY)'
            )
            time_values = [value.zfill(2) if index > 0 else value for index, value in enumerate(time_values)]
            sundata = {
                'date': reduce(lambda x, y: x + ':' + y, time_values),
                'lat': hdf['site_info']['site_lat'],
                'lon': hdf['site_info']['site_long'],
                'shadow_angle': (
                    -float(sdf[sample]['shadow_angle'])
                    if (
                        options.orientation_mode == 'block' and
                        options.reverse_block_sun_compass and
                        not _is_missing(sdf[sample]['shadow_angle'])
                    )
                    else sdf[sample]['shadow_angle']
                ),
                'delta_u': sdf[sample]['GMT_offset'],
            }
            sun_core_strike = sundec(sundata)
            if options.orientation_mode == 'block' and options.add_ninety_to_block_strike:
                sun_core_strike = _normalize_degrees(float(sun_core_strike) + 90.0)
            df[sample]['sun_core_strike'] = round(float(sun_core_strike), 1)

        igrf_required_fields = ['GMT_offset', 'year', 'month', 'days']
        if any(_is_missing(sdf[sample][field]) for field in igrf_required_fields):
            raise ValueError(
                'not enough data to calculate IGRF to correct bedding; '
                'input at least GMT_offset, year, month, and day of measurement'
            )

        if _is_missing(hdf['site_info']['site_elevation']):
            hdf['site_info']['site_elevation'] = 0.0

        for time_type in time_types:
            if _is_missing(sdf[sample][time_type]):
                sdf[sample][time_type] = 1

        date = to_year_fraction(
            dt(
                int(sdf[sample]['year']),
                int(sdf[sample]['month']),
                int(sdf[sample]['days']),
                int(sdf[sample]['hours']),
                int(sdf[sample]['minutes']),
            )
        )
        df[sample]['calculated_IGRF'] = list(
            igrf(
                [
                    date,
                    float(hdf['site_info']['site_elevation']) / 1000,
                    float(hdf['site_info']['site_lat']),
                    float(hdf['site_info']['site_long']),
                ]
            )
        )
        if float(df[sample]['calculated_IGRF'][0]) > 180:
            df[sample]['IGRF_local_dec'] = df[sample]['calculated_IGRF'][0] - 360
        else:
            df[sample]['IGRF_local_dec'] = df[sample]['calculated_IGRF'][0]

        logger(f"{hdf['site_info']['site_id']}{sample} has local IGRF declination of: ")
        logger(str(df[sample]['IGRF_local_dec']))

        logger('The local declination calculated through magnetic and sun compass comparison is:')
        if _is_missing(df[sample]['sun_core_strike']) or _is_missing(df[sample]['magnetic_core_strike']):
            df[sample]['calculated_mag_dec'] = 'insufficient data'
            logger('insufficient data')
        else:
            calc_mag_dec = (
                float(df[sample]['sun_core_strike']) -
                float(df[sample]['magnetic_core_strike'])
            )
            if calc_mag_dec > 180:
                df[sample]['calculated_mag_dec'] = calc_mag_dec - 360
            else:
                df[sample]['calculated_mag_dec'] = calc_mag_dec
            logger(f"    {df[sample]['calculated_mag_dec']:+.2f}")
            if abs(float(df[sample]['IGRF_local_dec']) - float(df[sample]['calculated_mag_dec'])) > 5:
                logger(
                    'WARNING: local IGRF declination & calculated magnetic declination '
                    'are more than 5 degree different'
                )
        logger('')

    logger('')
    logger('Site averages:')
    logger('Average of local IGRF declination is: ' + str(df.transpose()['IGRF_local_dec'].mean()))
    logger('')
    logger('---------------------OUTPUT-----------------------')


def _apply_orientation_preferences(
    hdf: pd.DataFrame,
    df: pd.DataFrame,
    options: CorrectionOptions,
) -> None:
    samples = df.keys()
    magnetic_declination_enabled = options.magnetic_declination_enabled()
    orientation_label = options.orientation_label()

    for sample in samples:
        if _is_missing(df[sample]['correct_bedding_using_local_dec']):
            df[sample]['correct_bedding_using_local_dec'] = 'yes'

        igrf_local_dec = _to_float(df[sample]['IGRF_local_dec'])
        sun_core_strike = _to_float(df[sample]['sun_core_strike'])
        magnetic_core_strike = _to_float(df[sample]['magnetic_core_strike'])

        if math.isnan(sun_core_strike) and math.isnan(magnetic_core_strike):
            raise ValueError(
                f'Sample {hdf["site_info"]["site_id"]}{sample} is missing both sun and magnetic orientation data.'
            )

        use_sun = options.prefer_sun_compass and not math.isnan(sun_core_strike)
        if use_sun:
            df[sample]['core_strike'] = round(float(sun_core_strike), 1)
            df[sample]['comment'] = 'sun compass orientation'
        else:
            if math.isnan(magnetic_core_strike):
                raise ValueError(
                    f'Sample {hdf["site_info"]["site_id"]}{sample} is missing magnetic orientation data.'
                )
            if magnetic_declination_enabled and not math.isnan(igrf_local_dec):
                df[sample]['core_strike'] = round(
                    _normalize_degrees(float(magnetic_core_strike) + float(igrf_local_dec)),
                    1,
                )
                df[sample]['comment'] = f'mag compass orientation ({orientation_label}, IGRF corrected)'
            else:
                df[sample]['core_strike'] = round(_normalize_degrees(float(magnetic_core_strike)), 1)
                df[sample]['comment'] = f'mag compass orientation ({orientation_label}, uncorrected)'

        should_correct_bedding = (
            options.apply_declination_to_bedding_strike and
            _yes(df[sample]['correct_bedding_using_local_dec'], default=True) and
            not _is_missing(df[sample]['bedding_strike']) and
            not math.isnan(igrf_local_dec)
        )
        if should_correct_bedding:
            df[sample]['corrected_bedding_strike'] = round(
                _normalize_degrees(float(df[sample]['bedding_strike']) + float(igrf_local_dec)),
                1,
            )


def _build_sam_header(hdf: pd.DataFrame, df: pd.DataFrame) -> str:
    site_values = ['site_lat', 'site_long']
    sam_header = hdf['site_info']['site_name'] + WINDOWS_NEWLINE

    for value in site_values:
        hdf['site_info'][value] = str(round(float(hdf['site_info'][value]), 1))
        if value == 'site_lat':
            sam_header += ' ' + hdf['site_info'][value]
        if value == 'site_long':
            sam_header += ' {:05.1f}'.format(float(hdf['site_info'][value]) % 360)
    sam_header += ' ' * 3 + '0.0'
    sam_header += WINDOWS_NEWLINE

    for sample in df.keys():
        sam_header += hdf['site_info']['site_id'] + str(sample) + WINDOWS_NEWLINE

    return sam_header


def _build_sample_file(sample: object, site_id: str, df: pd.DataFrame, options: CorrectionOptions, logger: Logger) -> str:
    attributes = ['core_strike', 'core_dip', 'bedding_strike', 'bedding_dip', 'mass']
    runs_value = df[sample]['runs']
    if _is_missing(runs_value):
        runs = []
    else:
        runs = str(runs_value).split(';')

    comment = df[sample]['comment']
    if _is_missing(comment):
        comment = ''

    assert len(site_id) <= 5, (
        'Locality ID exceeds 5 characters: refer to:'
        'http://cires.colorado.edu/people/jones.craig/PMag_Formats.html '
        '(although that says that 4 is the limit)'
    )
    assert len(comment) <= 255, (
        'Sample comment exceeds 255 characters: refer to:'
        'http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    )
    assert len(str(sample)) <= 9, (
        'Sample name exceeds 9 characters: refer to:'
        'http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    )

    new_file = site_id + ' ' + str(sample) + ' ' + str(comment) + WINDOWS_NEWLINE

    strat_level = df[sample]['strat_level']
    if _is_missing(strat_level):
        strat_level = '     0'
    strat_level = str(strat_level)
    assert len(strat_level) <= 6, (
        'Length of strat_level exceeds 6 characters: refer to:'
        'http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
    )
    new_file += ' ' + ' ' * (6 - len(strat_level)) + strat_level

    bedding_correction_enabled = (
        options.apply_declination_to_bedding_strike and
        _yes(df[sample]['correct_bedding_using_local_dec'], default=True) and
        not _is_missing(df[sample]['corrected_bedding_strike'])
    )

    for attribute in attributes:
        attribute_name = attribute
        attribute_value = df[sample][attribute_name]

        if attribute_name == 'bedding_strike' and _is_missing(attribute_value):
            attribute_value = 90.0
            df[sample][attribute_name] = attribute_value
        if attribute_name == 'bedding_dip' and _is_missing(attribute_value):
            attribute_value = 0.0
            df[sample][attribute_name] = attribute_value
        if attribute_name == 'bedding_strike' and bedding_correction_enabled:
            attribute_name = 'corrected_bedding_strike'
            attribute_value = df[sample][attribute_name]

        if _is_missing(attribute_value):
            if attribute_name == 'mass':
                attribute_string = '1.0'
                logger(f'no mass found for sample {sample}, setting to default = 1.0 g')
            else:
                attribute_string = ''
        else:
            attribute_string = str(round(float(attribute_value), 1))

        assert len(attribute_string) <= 5, (
            'Length of ' + attribute_name + ' exceeds 5 characters: refer to:' +
            'http://cires.colorado.edu/people/jones.craig/PMag_Formats.html'
        )

        new_file += ' ' + ' ' * (5 - len(attribute_string)) + attribute_string

    new_file += WINDOWS_NEWLINE

    for run in runs:
        if run:
            new_file += run + WINDOWS_NEWLINE

    return new_file.rstrip('\r\n') + WINDOWS_NEWLINE


def _field_magic_code(comments: pd.Series) -> str:
    normalized_comments = [str(comment).lower() for comment in comments.fillna('')]
    all_sun = all('sun compass orientation' in comment for comment in normalized_comments)
    all_mag = all('mag compass orientation' in comment for comment in normalized_comments)
    if all_sun:
        return 'SO-SUN'
    if all_mag:
        return 'SO-MAG'
    return 'SO-SM'


def generate_inp_file(output_directory: Path, df: pd.DataFrame, hdf: pd.DataFrame) -> Path:
    inps = ''
    inps += 'CIT' + WINDOWS_NEWLINE
    inps += (
        'sam_path\tfield_magic_codes\tlocation\tnaming_convention\t'
        'num_terminal_char\tdont_average_replicate_measurements\tpeak_AF\ttime_stamp' +
        WINDOWS_NEWLINE
    )
    inps += os.path.join('.', hdf['site_info']['site_id'] + '.sam') + '\t'
    inps += _field_magic_code(df.T['comment']) + '\t'

    site_name = hdf['site_info']['site_name']
    if _is_missing(site_name) or site_name == '':
        inps += 'unknown\t'
    else:
        inps += str(site_name) + '\t'

    first_sample_id = str(df.keys()[0])
    if first_sample_id[0] == '-' or hdf['site_info']['site_id'][-1] == '-':
        inps += '2\t'
        if first_sample_id[0] == '-':
            first_sample_id = first_sample_id.replace('-', '', 1)
    elif first_sample_id[0] == '.' or hdf['site_info']['site_id'][-1] == '.':
        inps += '3\t'
        if first_sample_id[0] == '.':
            first_sample_id = first_sample_id.replace('.', '', 1)
    elif hdf['site_info']['site_id'] == first_sample_id:
        inps += '5\t'
    else:
        inps += '4\t'

    sample_list = list(map(str, df.keys()))
    sample_ct = len(sample_list)
    char_num = len(min(sample_list, key=len))
    term_ct, term_unique = 0, 0
    the_rest = sample_list
    while term_ct < char_num:
        lastchar = [sample_name[-1] for sample_name in the_rest]
        the_rest = [sample_name[0:-1] for sample_name in the_rest]
        term_ct += 1
        unique_chars = np.unique(lastchar)
        unique_rest = np.unique(the_rest)
        if len(unique_chars) == sample_ct and len(unique_rest) == 1:
            term_unique = 1
            break
    inps += str(term_ct if term_unique else 0) + '\t'
    inps += '1\t0\t' + dt.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ') + WINDOWS_NEWLINE

    output_path = output_directory / (hdf['site_info']['site_id'] + '.inp')
    _write_text(output_path, inps)
    return output_path


def _write_updated_csv(
    file_name: str,
    output_directory: Path,
    hdf: pd.DataFrame,
    df: pd.DataFrame,
    sdf: pd.DataFrame,
    logger: Logger,
    output_name: str | None = None,
) -> Path:
    csv_str = ''
    csv_contents = _read_text_with_fallback(file_name).replace('\r\n', '\n').replace('\r', '\n')
    csv_lines = csv_contents.split('\n')

    for index in range(5):
        csv_str += csv_lines[index].rstrip('\n') + WINDOWS_NEWLINE

    comma_count = csv_lines[5].count(',')
    csv_str += 'site_elevation,' + str(hdf['site_info']['site_elevation']) + ',' * (comma_count - 1) + WINDOWS_NEWLINE

    header_line = csv_lines[6].rstrip('\n')
    csv_str += header_line + WINDOWS_NEWLINE
    header = header_line.split(',')
    samples = list(df.keys())

    for sample_index, sample in enumerate(samples):
        items = csv_lines[7 + sample_index].split(',')
        if len(items) < len(header):
            items.extend([''] * (len(header) - len(items)))
        items = [item.rstrip('\n') for item in items]
        for item_index in range(1, len(header)):
            column_name = header[item_index]
            if column_name == 'calculated_IGRF':
                value = df[sample][column_name]
                if isinstance(value, list):
                    items[item_index] = str(list(value)).replace(',', ';')
                else:
                    items[item_index] = str(value)
            elif column_name in df[sample].keys():
                items[item_index] = str(df[sample][column_name])
            elif column_name in sdf[sample].keys():
                items[item_index] = str(sdf[sample][column_name])
            else:
                raise KeyError('there is no item: ' + column_name)
        csv_str += reduce(lambda x, y: x + ',' + y, items) + WINDOWS_NEWLINE

    output_path = output_directory / (output_name or (hdf['site_info']['site_id'] + '.csv'))
    logger('Writing file - ' + str(output_path))
    _write_text(output_path, csv_str)
    return output_path


def _build_block_stitched_rows(source_path: Path, hdf: pd.DataFrame, df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    site_id = str(hdf['site_info']['site_id'])
    for sample in df.keys():
        rows.append(
            {
                'source_file': source_path.name,
                'site_id': site_id,
                'sample_name': str(sample),
                'comment': _stringify_csv_value(df[sample]['comment']),
                'magnetic_core_strike': _stringify_csv_value(df[sample]['magnetic_core_strike']),
                'sun_core_strike': _stringify_csv_value(df[sample]['sun_core_strike']),
                'IGRF_local_dec': _stringify_csv_value(df[sample]['IGRF_local_dec']),
                'calculated_mag_dec': _stringify_csv_value(df[sample]['calculated_mag_dec']),
                'corrected_block_strike': _stringify_csv_value(df[sample]['core_strike']),
                'core_dip': _stringify_csv_value(df[sample]['core_dip']),
                'bedding_strike': _stringify_csv_value(df[sample]['bedding_strike']),
                'bedding_dip': _stringify_csv_value(df[sample]['bedding_dip']),
            }
        )
    return rows


def write_stitched_block_orientation_file(
    output_directory: Path,
    rows: list[dict[str, object]],
    logger: Logger | None = None,
    file_name: str = 'all_blocks_corrected_orientations.csv',
) -> Path:
    logger = logger or _default_logger
    output_directory.mkdir(parents=True, exist_ok=True)
    output_path = output_directory / file_name

    frame = pd.DataFrame(rows)
    csv_text = frame.to_csv(index=False, lineterminator=WINDOWS_NEWLINE) if not frame.empty else ''
    logger('Writing file - ' + str(output_path))
    _write_text(output_path, csv_text)
    return output_path


def convert_csv_to_block_outputs(
    file_name: str,
    output_directory: str | None = None,
    options: CorrectionOptions | None = None,
    logger: Logger | None = None,
) -> ConversionResult:
    options = options or CorrectionOptions(orientation_mode='block')
    logger = logger or _default_logger

    source_path = Path(file_name).resolve()
    target_directory = Path(output_directory).resolve() if output_directory else source_path.parent
    target_directory.mkdir(parents=True, exist_ok=True)

    logger('Reading in file - ' + str(source_path))
    fix_line_breaks(str(source_path))
    hdf, df, sdf = _load_input_frames(str(source_path))
    _calculate_declinations(hdf, df, sdf, options, logger)
    _apply_orientation_preferences(hdf, df, options)

    site_id = str(hdf['site_info']['site_id'])
    corrected_csv = _write_updated_csv(
        str(source_path),
        target_directory,
        hdf,
        df,
        sdf,
        logger,
        output_name=f'{site_id}_block_corrected.csv',
    )
    stitched_rows = _build_block_stitched_rows(source_path, hdf, df)
    return ConversionResult(
        output_directory=target_directory,
        site_id=site_id,
        generated_files=[corrected_csv],
        stitched_rows=stitched_rows,
    )


def convert_csv_to_sam(
    file_name: str,
    output_directory: str | None = None,
    options: CorrectionOptions | None = None,
    logger: Logger | None = None,
) -> ConversionResult:
    options = options or CorrectionOptions()
    logger = logger or _default_logger

    source_path = Path(file_name).resolve()
    target_directory = Path(output_directory).resolve() if output_directory else source_path.parent
    target_directory.mkdir(parents=True, exist_ok=True)

    logger('Reading in file - ' + str(source_path))
    fix_line_breaks(str(source_path))
    hdf, df, sdf = _load_input_frames(str(source_path))
    _calculate_declinations(hdf, df, sdf, options, logger)
    _apply_orientation_preferences(hdf, df, options)

    site_id = str(hdf['site_info']['site_id'])
    sam_path = target_directory / f'{site_id}.sam'
    logger('Writing file - ' + str(sam_path))
    _write_text(sam_path, _build_sam_header(hdf, df))

    generated_files = [sam_path]
    for sample in df.keys():
        sample_path = target_directory / f'{site_id}{sample}'
        logger('Writing file - ' + str(sample_path))
        _write_text(sample_path, _build_sample_file(sample, site_id, df, options, logger))
        generated_files.append(sample_path)

    generated_files.append(_write_updated_csv(str(source_path), target_directory, hdf, df, sdf, logger))
    generated_files.append(generate_inp_file(target_directory, df, hdf))
    return ConversionResult(output_directory=target_directory, site_id=site_id, generated_files=generated_files)


def cli_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Create .sam headers and sample files from a site CSV template.')
    parser.add_argument('file_name', help='Input CSV file')
    parser.add_argument('output_directory', nargs='?', default=None, help='Optional output directory')
    parser.add_argument('--orientation-mode', choices=['core', 'block'], default='core')
    parser.add_argument('--prefer-sun-compass', dest='prefer_sun_compass', action='store_true', default=True)
    parser.add_argument('--ignore-sun-compass', dest='prefer_sun_compass', action='store_false')
    parser.add_argument('--correct-core-strike', dest='correct_core_strike', action='store_true', default=True)
    parser.add_argument('--no-correct-core-strike', dest='correct_core_strike', action='store_false')
    parser.add_argument('--correct-block-strike', dest='correct_block_strike', action='store_true', default=False)
    parser.add_argument('--no-correct-block-strike', dest='correct_block_strike', action='store_false')
    parser.add_argument('--correct-bedding-strike', dest='correct_bedding_strike', action='store_true', default=True)
    parser.add_argument('--no-correct-bedding-strike', dest='correct_bedding_strike', action='store_false')
    args = parser.parse_args(argv)

    options = CorrectionOptions(
        orientation_mode=args.orientation_mode,
        apply_declination_to_core_strike=args.correct_core_strike,
        apply_declination_to_block_strike=args.correct_block_strike,
        apply_declination_to_bedding_strike=args.correct_bedding_strike,
        prefer_sun_compass=args.prefer_sun_compass,
    )

    if args.orientation_mode == 'block':
        convert_csv_to_block_outputs(
            args.file_name,
            output_directory=args.output_directory,
            options=options,
        )
    else:
        convert_csv_to_sam(
            args.file_name,
            output_directory=args.output_directory,
            options=options,
        )
    return 0


__all__ = [
    'ConversionResult',
    'CorrectionOptions',
    'WINDOWS_NEWLINE',
    'cli_main',
    'convert_csv_to_block_outputs',
    'convert_csv_to_sam',
    'fix_line_breaks',
    'generate_inp_file',
    'write_stitched_block_orientation_file',
]