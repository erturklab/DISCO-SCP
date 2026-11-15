# Adapted from CosciaLab/Qupath_to_LMD (https://github.com/CosciaLab/Qupath_to_LMD, GPL-3.0).
# Changes: 8-well strip layout instead of a 384-well plate; Streamlit dependency removed (messages are printed, errors raised);
# shapes are written in the order of the samples-and-wells mapping; the calibration-triangle check was dropped.
# (See ../README.md for the full list of changes.)
import geopandas as gpd
import shapely
from lmd.lib import Collection
import ast
import pandas as pd
import numpy as np
from pathlib import Path
from loguru import logger

from qupath_to_lmd.utils import create_list_of_acceptable_wells, sample_placement_custom_strip, extract_coordinates


def load_and_QC_geojson_file(geojson_path, list_of_calibpoint_names=['calib1','calib2','calib3']):
    """
    Loads a GeoJSON file and checks for common issues before LMD conversion.
    
    Parameters
    ----------
    geojson_path : str or Path
        Path to the GeoJSON file
    list_of_calibpoint_names : list
        List of calibration point names
    
    Returns
    -------
    df : GeoDataFrame
        Filtered GeoDataFrame ready for processing
    calib_np_array : np.ndarray
        Numpy array of calibration points
    """

    df = gpd.read_file(geojson_path)
    logger.info(f"GeoJSON file loaded with shape {df.shape}")
    
    # Ensure 'name' column exists
    if 'name' not in df.columns:
        raise ValueError("No 'name' column found in GeoJSON; at least calibration points should be named")
    
    # Show geometry counts
    geometry_counts = df.geometry.geom_type.value_counts()
    for geom_type, count in geometry_counts.items():
        print(f"{geom_type}: {count}")
    
    # Check calibration points exist
    for point_name in list_of_calibpoint_names:
        if point_name not in df['name'].unique():
            print(f"WARNING: Calibration point '{point_name}' not found in GeoJSON")
    
    # Create calibration points array
    calib_np_array = np.array([
        [df.loc[df['name']==point_name, 'geometry'].values[0].x,
         df.loc[df['name']==point_name, 'geometry'].values[0].y]
        for point_name in list_of_calibpoint_names
    ])
    print(f"Calibration points array:\n{calib_np_array}")
    
    # Remove Point geometries
    df = df[df['geometry'].apply(lambda geom: not isinstance(geom, shapely.geometry.Point))]
    print(f"Removed point geometries. Remaining shapes: {len(df)}")
    
    # Remove empty classifications
    n_nan = df['classification'].isna().sum()
    if n_nan > 0:
        print(f"{n_nan} unclassified objects will be ignored")
        df = df[df['classification'].notna()]
    
    # Extract classification name
    df['classification_name'] = df['classification'].apply(
        lambda x: ast.literal_eval(x).get('name') if isinstance(x, str) else x.get('name')
    )
    
    # Remove MultiPolygons
    n_multipoly = (df.geometry.geom_type == 'MultiPolygon').sum()
    if n_multipoly > 0:
        print(f"{n_multipoly} MultiPolygon objects will be ignored")
        df = df[df.geometry.geom_type != 'MultiPolygon']
    
    # Simplify coordinates
    df['coords'] = df.geometry.simplify(1).apply(lambda g: np.array(g.exterior.coords) if isinstance(g, shapely.geometry.Polygon) else None)
    
    print("GeoJSON QC complete.")
    return df, calib_np_array


def load_and_QC_SamplesandWells(geojson_path, list_of_calibpoint_names, samples_and_wells_input):
    """
    Load a GeoJSON and a samples-and-wells mapping, check for common issues.

    Parameters
    ----------
    geojson_path : str or Path
        Path to the GeoJSON file.
    list_of_calibpoint_names : list
        List of calibration point names to exclude from samples.
    samples_and_wells_input : dict
        Mapping of sample names to wells.

    Returns
    -------
    samples_and_wells : dict
        The validated samples and wells mapping.
    """

    # Load GeoJSON
    df = gpd.read_file(geojson_path)

    # Remove calibration points
    df = df[~df['name'].isin(list_of_calibpoint_names)]

    # Keep only objects with classification and exclude MultiPolygons
    df = df[df['classification'].notna()]
    df = df[df.geometry.geom_type != 'MultiPolygon']

    # Extract sample names from classification
    df['Name'] = df['classification'].apply(
        lambda x: ast.literal_eval(x).get('name') if isinstance(x, str) else x.get('name')
    )

    # Use input dict directly
    if not isinstance(samples_and_wells_input, dict):
        raise ValueError("samples_and_wells_input must be a dict mapping sample names to wells")
    samples_and_wells = samples_and_wells_input

    # Check wells against acceptable wells
    list_of_acceptable_wells = create_list_of_acceptable_wells()
    warning_limit = 0
    for well in samples_and_wells.values():
        if well not in list_of_acceptable_wells:
            warning_limit += 1
            if warning_limit <= 10:
                print(f"WARNING: Well '{well}' is not in the list of acceptable wells (rows A-H, columns 1-24)")

    if warning_limit > 10:
        print("WARNING: More than 10 wells are not acceptable. Check your mapping!")

    # Check sample names against GeoJSON
    unmatched_names = [name for name in df['Name'].unique() if name not in samples_and_wells.keys()]
    if unmatched_names:
        print("WARNING: These sample names are in the GeoJSON but not in your mapping:")
        print(unmatched_names)

    print("Samples and wells QC done!")

    return samples_and_wells

def create_collection(
    geojson_path: str = None,
    list_of_calibpoint_names: list = None,
    samples_and_wells_input: dict = None
):

    # Load GeoJSON
    df = gpd.read_file(geojson_path)

    # Extract calibration points
    calib_np_array = np.array([
        [
            df.loc[df['name'] == point_name, 'geometry'].values[0].x,
            df.loc[df['name'] == point_name, 'geometry'].values[0].y
        ]
        for point_name in list_of_calibpoint_names
    ])

    # Filter GeoJSON: exclude points, null classifications, multipolygons
    df = df[df['geometry'].apply(lambda geom: not isinstance(geom, shapely.geometry.Point))]
    df = df[df['classification'].notna()]
    df = df[df.geometry.geom_type != 'MultiPolygon']

    # Extract coordinates and sample names
    df['coords'] = df.geometry.simplify(1).apply(extract_coordinates)
    df['Name'] = df['classification'].apply(
        lambda x: ast.literal_eval(x).get('name') if isinstance(x, str) else x.get('name')
    )

    # Validate input
    if samples_and_wells_input is None or not isinstance(samples_and_wells_input, dict):
        raise ValueError("Please provide a samples_and_wells_input dictionary")
    samples_and_wells = samples_and_wells_input

    # Filter out shapes not in mapping
    df = df[df['Name'].isin(samples_and_wells.keys())]

    # --- SORT according to mapping order to match CapIDs ---
    mapping_order = list(samples_and_wells.keys())
    df['mapping_index'] = df['Name'].apply(lambda x: mapping_order.index(x))
    df = df.sort_values('mapping_index').reset_index(drop=True)

    # Create LMD collection
    the_collection = Collection(calibration_points=calib_np_array)
    the_collection.orientation_transform = np.array([[1, 0], [0, -1]])

    for i in df.index:
        shape_name = df.at[i, 'Name']
        well = samples_and_wells[shape_name]
        the_collection.new_shape(df.at[i, 'coords'], well=well)

    # Plot and save
    the_collection.plot(save_name="./TheCollection.png")
    print("Collection plotted: ./TheCollection.png")
    
    try:
        stats = the_collection.stats()
        print(stats)
    except Exception as e:
        logger.warning(f"Could not calculate stats: {e}")

    # Save XML
    the_collection.save(f'./{Path(geojson_path).name.replace("geojson", "xml")}')
    print(f"Collection saved: ./{Path(geojson_path).name.replace('geojson', 'xml')}")

    # Export 8-well strip CSV
    df_strip = sample_placement_custom_strip(samples_and_wells)
    csv_name = f'./{Path(geojson_path).name.replace("geojson", "_8well_strip.csv")}'
    df_strip.to_csv(csv_name, index=False)
    print(f"Strip CSV saved: {csv_name}")

    return the_collection, df_strip