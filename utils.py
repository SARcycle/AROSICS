def get_credentials(tag):
    import json
    with open('/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/credentials.json') as f:
        credentials = json.load(f)
        try:
            username = credentials[tag]['username']
        except:
            username = None
        try:
            password = credentials[tag]['password']
        except:
            password = None

    return username, password


def extract_specific_files(zip_path, target_folder, search_str, target_filename):
    import zipfile, os, re
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for file in zip_ref.infolist():
            if re.search(search_str, file.orig_filename):
                # Extract the file to a temporary path
                file.filename = os.path.basename(file.filename)
                temp_path = os.path.join(target_folder, file.filename)
                zip_ref.extract(file, target_folder)
                # Move the file to the target directory, without subdirectories
                os.rename(temp_path, os.path.join(target_folder, target_filename))


def get_epsg(tif_file):
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
    import xml.etree.ElementTree as ET

    # Load the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()
    for elem in root.iter():
        if 'name' in elem.attrib and elem.attrib['name'] == key:
            return float(elem.text)
    return None
