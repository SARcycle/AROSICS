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
