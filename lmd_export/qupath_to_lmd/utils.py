# Adapted from CosciaLab/Qupath_to_LMD (https://github.com/CosciaLab/Qupath_to_LMD, GPL-3.0).
# Changes: well placement for an 8-well strip instead of a 384-well plate; Streamlit dependency and unused functions removed.
# (See ../README.md for the full list of changes.)
import string
import json
from pathlib import Path

import numpy as np
import pandas as pd


def create_list_of_acceptable_wells():
    """Return all wells in column-major order for an 8-row (A-H) plate with 24 columns."""
    rows = list(string.ascii_uppercase[:8])
    cols = range(1, 25)
    return [f"{row}{col}" for col in cols for row in rows]


def sample_placement_custom_strip(samples_and_wells):
    """Return a table (well, sample, row, col) for the sample-to-well mapping."""
    rows, cols, sample_names = [], [], []
    for sample, well in samples_and_wells.items():
        rows.append(well[0])
        cols.append(int(well[1:]))
        sample_names.append(sample)
    return pd.DataFrame({'well': list(samples_and_wells.values()),
                         'sample': sample_names,
                         'row': rows,
                         'col': cols})


def create_default_wells_column_major(num_wells):
    """Generate wells in column-major order for 8 rows (A-H) up to num_wells."""
    rows = list(string.ascii_uppercase[:8])
    wells = []
    col = 1
    while len(wells) < num_wells:
        for row in rows:
            wells.append(f"{row}{col}")
            if len(wells) >= num_wells:
                break
        col += 1
    return wells


def assign_samples_to_wells(categories1, categories2, n_reps,
                            output_folder=".", filename="samples_and_wells.txt"):
    """Create the sample-to-well mapping (column-major: A1-H1, A2-H2, ...) and save it as JSON text."""
    num_samples = len(categories1) * len(categories2) * n_reps
    wells = create_default_wells_column_major(num_samples)
    mapping = {}

    sample_idx = 0
    for cat1 in categories1:
        for cat2 in categories2:
            for rep in range(1, n_reps + 1):
                mapping[f"{cat1}_{cat2}_{rep}"] = wells[sample_idx]
                sample_idx += 1

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    output_path = output_folder / filename
    output_path.write_text(json.dumps(mapping, indent=4))

    print(f"Saved mapping to {output_path}")
    return mapping


def extract_coordinates(geom):
    """Return a numpy array of coordinates from a shapely Polygon or LineString."""
    if geom.geom_type == "Polygon":
        return np.array(geom.exterior.coords)
    elif geom.geom_type == "LineString":
        return np.array(geom.coords)
    else:
        raise ValueError(f"Unsupported geometry type: {geom.geom_type}")
