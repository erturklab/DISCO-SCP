# QuPath annotations to Leica LMD XML

Converts QuPath annotations (GeoJSON) into the cutting file (XML) for the Leica LMD, with the shapes assigned to the wells of an 8-well strip. Used to collect the single cells and 10-cell areas of the study.

This folder is adapted from [CosciaLab/Qupath_to_LMD](https://github.com/CosciaLab/Qupath_to_LMD) (GPL-3.0), which was written for a 384-well plate. The XML itself is generated with [py-lmd](https://github.com/MannLabs/py-lmd) (Apache 2.0).

**Changes to the original** (compared with the Qupath_to_LMD version of September 2025, which is the closest match to the code we started from):

1. **Wells:** shapes are placed in an 8-well strip (column-major, A1–H1, A2–H2, ...) instead of a 384-well plate, and the well table is written as `<name>._8well_strip.csv`.
2. **Order of the shapes:** the shapes are written to the XML in the order of the samples-and-wells list (sample 1, sample 2, ...), so that the cutting order matches the order of the collection wells. The original writes them in the order of the GeoJSON file.
3. **No web app:** the Streamlit interface was removed. Checks print messages in the notebook, and errors stop the notebook with an error message.
4. **Dropped functions:** the CSV metadata colouring, the class-name generator and the 384-well plate functions are not included.
5. **Dropped check:** the original warned if fewer than 25% of the shapes lie inside the triangle of the three calibration points. This check is not included here, so check the overview image (`TheCollection.png`) by eye.

## Contents

| File | What it does |
|---|---|
| `notebooks/geojson_to_lmd_xml.ipynb` | Step-by-step notebook: check the GeoJSON, assign wells, create the XML |
| `qupath_to_lmd/utils.py` | Well assignment and well table |
| `qupath_to_lmd/qc_utils.py` | Checks of the GeoJSON and of the shape-to-well mapping, creation of the XML (`create_collection`) |
| `requirements.txt` | Python packages |

## Input from QuPath

1. Draw the annotations and give each one a class name (e.g. `Somatosensory_Cortex_Cell_1`).
2. Add three calibration points with the point tool, named for example `calib_1`, `calib_2`, `calib_3`.
3. Export the annotations as a GeoJSON FeatureCollection.

## How to run

1. Create a Python environment and install the packages: `pip install -r requirements.txt`
2. Open `notebooks/geojson_to_lmd_xml.ipynb`, set the paths and names in section 1, and run all cells.
3. Output in the folder you set: `<name>.xml` (load into the Leica LMD software), `TheCollection.png` (overview of the shapes), `<name>._8well_strip.csv` (well layout) and `Samples_and_wells.txt`.

The number of shapes in the GeoJSON must equal `n_reps` times the number of category combinations, otherwise the shape-to-well check lists the names that do not match.

## License and credit

GPL-3.0 (see `LICENSE`), because this folder is derived from GPL-3.0 code. If you use it, please also cite the original project: Nimo et al., OpenDVP, bioRxiv 2025, doi:10.1101/2025.07.13.662099, and py-lmd (MannLabs).
