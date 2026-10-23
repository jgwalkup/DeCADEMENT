"""
branch_output.py
================
Split output of branched runs (run_multi_climate_simulations.py with
branch_pulse > 0), and reading it back as full series.

Layout written by the runner (csv output, default for branched runs):

    <output-dir>/spinup/<community>_s<seed>/          days before the branch, once
        <Table>.csv                                   day columns 0 .. B (B = branch day)
        Initialization/                               complete
    <output-dir>/<scenario>/<community>_s<seed>/      one per climate scenario
        <Table>.csv                                   day columns B+1 .. pulse x end_time
        Initialization/Temp.csv, Psi.csv              the scenario's climate only
        branch_output.json                            what was split, and where the rest is

Day-indexed tables have one column per output day (labels 0 .. n_days, column
k = state after day k-1); Taxon_count_repop has one column per day with labels
restarting every pulse; Kill has one row per output day. Every other file is
written whole in each scenario folder.

For analysis, read tables through read_table(), which returns the full series
(spin-up + branch) for split runs and the file as is otherwise:

    import sys; sys.path.insert(0, 'analysis')
    from branch_output import read_table, initialization_dir
    microbes = read_table('sim_heatwave/climate_summer_warm/community_0001_s0', 'MicrobesSeries')
    init = initialization_dir('sim_heatwave/climate_summer_warm/community_0001_s0')
"""
import csv
import json
import os
from pathlib import Path

import pandas as pd

MARKER = 'branch_output.json'
ROW_DAY_TABLES = ('Kill',)          # tables with one row per output day
INIT_KEEP = ('Temp.csv', 'Psi.csv')  # Initialization files kept per scenario


def _is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def _day_columns(header, min_days):
    """
    (first day-column position, kind) for a day-indexed header, else None.
    kind 'daily': labels 0, 1, 2, ... (initial state + one column per day);
    kind 'per_pulse': labels restart each pulse (one column per day).
    """
    p = 0
    while p < len(header) and not _is_int(header[p]):
        p += 1
    labels = header[p:]
    if p == 0 or len(labels) < min_days or not all(_is_int(x) for x in labels):
        return None
    ints = [int(x) for x in labels]
    if ints == list(range(len(ints))):
        return p, 'daily'
    if ints[0] == 0 and any(b < a for a, b in zip(ints, ints[1:])):
        return p, 'per_pulse'
    return None


def _rewrite(path, rows):
    tmp = path.with_suffix('.csv.tmp')
    with open(tmp, 'w', newline='') as fh:
        csv.writer(fh, lineterminator='\n').writerows(rows)
    os.replace(tmp, path)


def split_run_output(run_dir, part, branch_day):
    """
    Trim a run folder written by Output.export_to_csv to one part of a branched
    run, in place, keeping the CSV text of every kept value unchanged.

    part 'spinup': day columns/rows before the branch (daily tables: initial
                   state + days 0 .. branch_day - 1).
    part 'branch': day columns/rows from the branch on; Initialization reduced
                   to Temp.csv and Psi.csv.
    Returns {'columns': [...], 'rows': [...]} naming the trimmed tables.
    """
    run_dir = Path(run_dir)
    trimmed = {'columns': [], 'rows': []}
    for f in sorted(run_dir.glob('*.csv')):
        with open(f, newline='') as fh:
            rows = list(csv.reader(fh))
        if not rows:
            continue
        name = f.stem
        if name in ROW_DAY_TABLES:
            body = rows[1:]
            cut = branch_day + 1          # rows: initial + days before the branch
            keep = body[:cut] if part == 'spinup' else body[cut:]
            _rewrite(f, [rows[0]] + keep)
            trimmed['rows'].append(name)
            continue
        found = _day_columns(rows[0], min_days=branch_day)
        if found is None:
            continue
        p, kind = found
        cut = p + (branch_day + 1 if kind == 'daily' else branch_day)
        if part == 'spinup':
            out = [r[:cut] for r in rows]
        else:
            out = [r[:p] + r[cut:] for r in rows]
        _rewrite(f, out)
        trimmed['columns'].append(name)

    if part == 'branch':
        init = run_dir / 'Initialization'
        if init.is_dir():
            for g in init.iterdir():
                if g.is_file() and g.name not in INIT_KEEP:
                    g.unlink()
    return trimmed


def write_marker(run_dir, spinup_dir, branch_day, n_days, trimmed):
    run_dir = Path(run_dir)
    info = {
        'spinup_output': os.path.relpath(spinup_dir, run_dir),
        'branch_day': int(branch_day),
        'n_days': int(n_days),
        'trimmed_columns': trimmed['columns'],
        'trimmed_rows': trimmed['rows'],
        'note': ('Day-indexed tables hold only days from the branch on; the '
                 'days before are in spinup_output (shared by all scenarios '
                 'of this community). Read full series with '
                 'branch_output.read_table().'),
    }
    with open(run_dir / MARKER, 'w') as fh:
        json.dump(info, fh, indent=1)


def _read_csv(source, header):
    n_index = 0
    while n_index < len(header) and header[n_index] == '':
        n_index += 1
    n_index = max(n_index, 1)
    return pd.read_csv(source, index_col=list(range(n_index)),
                       float_precision='round_trip')


def _rows(path):
    with open(path, newline='') as fh:
        return list(csv.reader(fh))


def branch_info(run_dir):
    """The branch_output.json contents of a split run folder, or None."""
    m = Path(run_dir) / MARKER
    if not m.exists():
        return None
    with open(m) as fh:
        return json.load(fh)


def read_table(run_dir, name):
    """
    One output table as a DataFrame: for split branched runs the spin-up days
    and the branch days joined into the full series; otherwise the file as is.
    """
    import io
    run_dir = Path(run_dir)
    path = run_dir / f'{name}.csv'
    info = branch_info(run_dir)
    split_cols = info is not None and name in info['trimmed_columns']
    split_rows = info is not None and name in info['trimmed_rows']
    if not (split_cols or split_rows):
        with open(path, newline='') as fh:
            header = next(csv.reader(fh))
        return _read_csv(path, header)
    # Rebuild the full CSV text from the two parts and parse it once, so the
    # result is exactly what reading an unsplit file gives (including pandas'
    # renaming of repeated column labels, e.g. in Taxon_count_repop).
    spin = (run_dir / info['spinup_output']).resolve()
    rs, rb = _rows(spin / f'{name}.csv'), _rows(path)
    if split_cols:
        p = _day_columns(rs[0], min_days=1)[0]
        rows = [x + y[p:] for x, y in zip(rs, rb)]
    else:
        rows = rs + rb[1:]
    buf = io.StringIO()
    csv.writer(buf, lineterminator='\n').writerows(rows)
    buf.seek(0)
    return _read_csv(buf, rows[0])


def initialization_dir(run_dir):
    """
    The folder holding the run's Initialization tables. For split runs this is
    the spin-up's (shared); the scenario's own Temp.csv and Psi.csv stay in
    <run_dir>/Initialization.
    """
    run_dir = Path(run_dir)
    info = branch_info(run_dir)
    if info is None:
        return run_dir / 'Initialization'
    return (run_dir / info['spinup_output']).resolve() / 'Initialization'
