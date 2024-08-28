import re
from datetime import datetime

import coreg_main


def get_credentials(tag):
    """
    Retrieves the username and password from the credentials file based on the provided tag.

    Args:
        tag (str): The tag to retrieve the credentials for.

    Returns:
        tuple: A tuple containing the username and password.
    """

    import json

    # Open the credentials file
    with open('/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/credentials.json') as f:
        # Load the credentials from the file
        credentials = json.load(f)

        # Initialize the username and password to None
        username = None
        password = None

        # Try to retrieve the username and password from the credentials
        try:
            username = credentials[tag]['username']
        except:
            pass
        try:
            password = credentials[tag]['password']
        except:
            pass

    # Return the username and password
    return username, password


def extract_specific_files(zip_path, target_folder, search_str, target_filename):
    """
    Extracts specific files from a zip archive based on a search string.

    Args:
        zip_path (str): The path to the zip archive.
        target_folder (str): The target folder to extract the files to.
        search_str (str): The search string to match files against.
        target_filename (str): The target filename to rename the extracted file to.

    Returns:
        None
    """

    import zipfile, os, re

    # Open the zip archive
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Iterate over the files in the zip archive
        for file in zip_ref.infolist():
            # Check if the file matches the search string
            if re.search(search_str, file.orig_filename):
                # Extract the file to a temporary path
                file.filename = os.path.basename(file.filename)
                temp_path = os.path.join(target_folder, file.filename)
                zip_ref.extract(file, target_folder)
                # Move the file to the target directory, without subdirectories
                os.rename(temp_path, os.path.join(target_folder, target_filename))


def get_epsg(tif_file):
    """
    Retrieves the EPSG code from a GeoTIFF file.

    Args:
        tif_file (str): The path to the GeoTIFF file.

    Returns:
        int: The EPSG code.
    """

    from osgeo import gdal, osr

    # Open the dataset
    dataset = gdal.Open(tif_file)
    # Get the spatial reference
    srs = osr.SpatialReference()
    srs.ImportFromWkt(dataset.GetProjection())
    # Get the EPSG code
    epsg_code = srs.GetAuthorityCode(None)
    return int(epsg_code)


def get_pixel_spacing(file_path):
    from osgeo import gdal
    import numpy as np

    # Open the raster file
    ds = gdal.Open(file_path)

    # Check if the file was opened successfully
    if ds is None:
        raise RuntimeError(f"Failed to open file: {file_path}")

    try:
        # Get the geotransform
        gt = ds.GetGeoTransform()

        # Get the pixel spacing (x and y directions)
        pixel_spacing_x = np.round(gt[1])
        pixel_spacing_y = np.round(-gt[5])

        return pixel_spacing_x, pixel_spacing_y
    finally:
        # Close the dataset to free up memory
        ds = None


def extract_value_from_xml(file_path, key):
    """
    Extracts a value from an XML file based on a key.

    Args:
        file_path (str): The path to the XML file.
        key (str): The key to extract the value for.

    Returns:
        float: The extracted value, or None if not found.
    """

    import xml.etree.ElementTree as ET

    # Load the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()
    for elem in root.iter():
        # Check if the element has a name attribute matching the key
        if 'name' in elem.attrib and elem.attrib['name'] == key:
            # Return the value as a float
            return float(elem.text)
    # Return None if the key is not found
    return None


def change_resolution_CRL(CRL, gsd_new=None, backup=None):
    """
    Changes the resolution within the CRL.tiepoint_grid.shift to a different GSD to minimize requirement in the interpolation step.
    CAUTION: since CRL.tiepoint_grid.shift is mutable, changes to CRL_tmp are also reflected in the original CRL -> only a reference, not a copy.


    Args:
        CRL: AROSICS COREG object.
        gsd (int/float): The desired gsd.
        backup: a backup of the changed parameters to be reinstated

    Returns:
        CRL: AROSICS COREG object.
        backup: a backup of the changed parameters to be saved
    """
    import copy
    import numpy as np

    class coreg_backup:
        def __init__(self, CRL):
            self.arr = copy.deepcopy(CRL.tiepoint_grid.shift.arr)
            self.geotransform = copy.deepcopy(CRL.tiepoint_grid.shift.geotransform)
            self.GCPList = []
            for gcp in CRL.tiepoint_grid.GCPList:
                self.GCPList.append([gcp.GCPLine, gcp.GCPPixel])
            self.CoRegPoints_table = copy.deepcopy(CRL.tiepoint_grid.CoRegPoints_table)

    if backup is None:
        backup = coreg_backup(CRL=CRL)
        img_dims = coreg_main.get_extent_and_dimensions(CRL.im2shift.filePath)
        gsd_old = CRL.tiepoint_grid.shift.geotransform[1]
        CRL.tiepoint_grid.shift.arr = np.zeros(
            (len(np.arange(img_dims[2], img_dims[3], gsd_new)), len(np.arange(img_dims[0], img_dims[1], gsd_new))),
            dtype=CRL.tiepoint_grid.shift.dtype)  # First y direction, then x direction
        CRL.tiepoint_grid.shift.geotransform[1] = gsd_new
        CRL.tiepoint_grid.shift.geotransform[5] = -1 * gsd_new
        for gcp in CRL.tiepoint_grid.GCPList:
            gcp.GCPPixel = np.round(gcp.GCPPixel / (gsd_new / gsd_old))
            gcp.GCPLine = np.round(gcp.GCPLine / (gsd_new / gsd_old))
        CRL.tiepoint_grid.CoRegPoints_table.X_IM = np.round(
            CRL.tiepoint_grid.CoRegPoints_table.X_IM / (gsd_new / gsd_old))
        CRL.tiepoint_grid.CoRegPoints_table.Y_IM = np.round(
            CRL.tiepoint_grid.CoRegPoints_table.Y_IM / (gsd_new / gsd_old))
        CRL.tiepoint_grid.CoRegPoints_table.X_SHIFT_PX = np.round(
            CRL.tiepoint_grid.CoRegPoints_table.X_SHIFT_PX / (gsd_new / gsd_old))
        CRL.tiepoint_grid.CoRegPoints_table.Y_SHIFT_PX = np.round(
            CRL.tiepoint_grid.CoRegPoints_table.Y_SHIFT_PX / (gsd_new / gsd_old))
        return CRL, backup

    else:
        CRL.tiepoint_grid.shift.arr = copy.deepcopy(backup.arr)
        CRL.tiepoint_grid.shift.geotransform = copy.deepcopy(backup.geotransform)
        for gcp_orig, gcp_backup in zip(CRL.tiepoint_grid.GCPList, backup.GCPList):
            gcp_orig.GCPLine = gcp_backup[0]
            gcp_orig.GCPPixel = gcp_backup[1]
        CRL.tiepoint_grid.CoRegPoints_table = backup.CoRegPoints_table
        return CRL


def determine_cores(tiles):
    import os
    # If there are four tiles or less use all cores
    if tiles <= 4:
        return os.cpu_count() - 1
    elif tiles >= 11:
        return min(24, os.cpu_count() - 1)
    else:
        return min(int((tiles - 4) / (11 - 4) * (24 - os.cpu_count()) + os.cpu_count()), os.cpu_count() - 1)


def parse_date(date_str):
    """
    Parse a date string in yyyy-mm-dd, yyyy/mm/dd, or dd.mm.yyyy format.

    Args:
        date_str (str): The date string to parse.

    Returns:
        datetime: The parsed date.
    """
    patterns = [
        r"(\d{4})-(\d{2})-(\d{2})",  # yyyy-mm-dd
        r"(\d{4})/(\d{2})/(\d{2})",  # yyyy/mm/dd
        r"(\d{2})\.(\d{2})\.(\d{4})",  # dd.mm.yyyy
        r"(\d{4})(\d{2})(\d{2})",  # yyyymmdd
    ]

    for pattern in patterns:
        match = re.match(pattern, date_str)
        if match:
            if len(match.groups()) == 3:
                year, month, day = map(int, match.groups())
                return datetime(year, month, day)

    raise ValueError("Invalid date format")
