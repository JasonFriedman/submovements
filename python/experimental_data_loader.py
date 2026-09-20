"""Dataset-specific loaders for experimental submovement data."""

from pathlib import Path
import re
import numpy as np
from scipy.signal import filtfilt, butter


def load_tablet_csv_directory(dir_name):
    """
    Load tablet experiment data from a directory of trial CSV files.

    Expected per-file format:
    - Columns 0 and 1: x/y position
    - Column 3: pen pressure (samples with pressure > 0 are retained)
    - Column 4: timestamp (milliseconds)

    Filenames are expected to match: ``tb_*block<block>_trial<trial>.csv``.

    Returns:
        position_filtered (list[np.ndarray]): Filtered position arrays, one per trial.
        velocity (list[np.ndarray]): Velocity arrays (same shape as position), one per trial.
        time (list[np.ndarray]): Time arrays in seconds, each starting at 0.
    """

    directory = Path(dir_name)
    csv_files = sorted(directory.glob('*.csv'))
    if not csv_files:
        raise ValueError('Must specify a directory to load the csv files from')

    trial_records = []
    for file_path in csv_files:
        match = re.search(r'tb_.*block(\d*)_trial(\d*)\.csv', file_path.name)
        if match is None:
            continue
        trial_records.append((int(match.group(1)), int(match.group(2)), file_path))

    if not trial_records:
        raise ValueError('No files matched expected pattern tb_*block*_trial*.csv')

    trial_records.sort(key=lambda item: (item[0], item[1]))

    position_filtered = []
    velocity = []
    time = []

    for _, _, file_path in trial_records:
        data = np.loadtxt(file_path, delimiter=',')
        pressure = data[:, 3]
        position = data[pressure > 0, :2] / 1000
        _time = data[pressure > 0, 4] / 1000

        if position.shape[0] < 3:
            continue

        _time = _time - _time[0]
        dt = np.median(np.diff(_time))
        b, a = butter(2, 5 / ((1 / dt) / 2))
        _position_filtered = filtfilt(b, a, position, axis=0)
        _velocity = np.vstack([[0, 0], np.diff(_position_filtered, axis=0) / dt])

        time.append(_time)
        position_filtered.append(_position_filtered)
        velocity.append(_velocity)

    return position_filtered, velocity, time
