import hashlib
import os.path
import sys
from datetime import datetime

import geopandas as gpd
import requests

import utils


def image_search(start_date, end_date):
    """
    Searches for satellite images within a specified date range and area of interest.
    Args:
        start_date (datetime): Start date for the search.
        end_date (datetime): End date for the search.
    Returns:
        dict: A dictionary containing search results.
            Keys: 'tile_0000', 'tile_0001', ...
            Values: {'download_id', 'acquisition_date', 'relative_orbit', 'scene_name', 'checksum'}
    """
    print(f"Searching for images")
    if isinstance(start_date, datetime):
        start_date = start_date.strftime("%Y-%m-%d")
    if isinstance(end_date, datetime):
        end_date = end_date.strftime("%Y-%m-%d")

    proc_lvl = "L2A"
    data_collection = "SENTINEL-2"

    # Read the geospatial data
    gdf = gpd.read_file(
        "/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/swissBOUNDARIES3D_1_5_WGS84_buffered_5000m_simplified_DP_5000m.gpkg",
        layer="bbox",
    )
    aoi = gdf["geometry"][0]

    # Construct the request URL
    request_url = (
        f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?"
        f"$filter=Collection/Name eq '{data_collection}' "
        f"and OData.CSC.Intersects(area=geography'SRID=4326;{aoi}') "
        f"and ContentDate/Start gt {start_date}T00:00:00.000Z "
        f"and ContentDate/Start lt {end_date}T00:00:00.000Z&$top=1000"
    )

    # Fetch data from the API
    json = requests.get(request_url).json()

    # Extract relevant information
    ids = [d["Id"] for d in json["value"] if "Id" in d and proc_lvl in d.get("Name", "")]
    acquisition_dates = [
        datetime.strptime(d["ContentDate"]["Start"], "%Y-%m-%dT%H:%M:%S.%fZ")
        for d in json["value"]
        if "ContentDate" in d and proc_lvl in d.get("Name", "")
    ]
    scene_names = [
        d["Name"] for d in json["value"] if "Name" in d and proc_lvl in d.get("Name", "")
    ]
    checksums = [
        d["Checksum"][0]["Value"]
        for d in json["value"]
        if "Checksum" in d and proc_lvl in d.get("Name", "")
    ]
    rel_orbits = [int(scene_name.split("_")[4].lstrip("R")) for scene_name in scene_names if
                  len(scene_name.split("_")) >= 5]

    # Create a dictionary for search results
    search_result = {}
    for i in range(len(ids)):
        search_result[f"tile_{i:04}"] = {
            "download_id": ids[i],
            "acquisition_date": acquisition_dates[i],
            "relative_orbit": rel_orbits[i],
            "scene_name": scene_names[i],
            "checksum": checksums[i],
        }

    return search_result

def get_access_token():
    """
    Retrieves an access token for Copernicus services.
    Returns:
        str: Access token.
    """
    username, password = utils.get_credentials("copernicus")
    data = {
        "client_id": "cdse-public",
        "username": username,
        "password": password,
        "grant_type": "password",
    }
    try:
        r = requests.post(
            "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
            data=data,
        )
        r.raise_for_status()
    except Exception as e:
        raise Exception(
            f"Access token creation failed. Response from the server was: {r.json()}"
        )
    return r.json()["access_token"]

def calculate_md5(file_path):
    """
    Calculates the MD5 hash of a file.
    Args:
        file_path (str): Path to the file.
    Returns:
        str: MD5 hash.
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def S2_scene_download(base_path, start_date=None, end_date=None, search_result=None):
    """
    Downloads Sentinel-2 scenes based on the provided start and end dates.

    Args:
        base_path (str): The base path to save the downloaded scenes.
        start_date (datetime, optional): The start date for the scene search. Defaults to None.
        end_date (datetime, optional): The end date for the scene search. Defaults to None.
        search_result (dict, optional): The search result containing the scene information. Defaults to None.
    """

    import requests
    from datetime import datetime
    from requests.packages.urllib3.exceptions import InsecureRequestWarning

    # If search_result is provided, use it
    if search_result is not None:
        pass
    # If start_date is not provided, exit
    elif start_date is None:
        print('No start date provided. Exiting.')
        sys.exit()
    # If end_date is not provided, set it to the current date
    elif end_date is None:
        end_date = datetime.now()
        print('No end date provided. Defaulting to current date.')
    else:
        # Disable insecure request warnings
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

        # Convert start_date and end_date to strings if they are datetime objects
        if isinstance(start_date, datetime):
            start_date = start_date.strftime('%Y-%m-%d')
        if isinstance(end_date, datetime):
            end_date = end_date.strftime('%Y-%m-%d')

        # Perform the scene search
        search_result = image_search(start_date=start_date, end_date=end_date)
        search_result = clean_search_result(search_result, base_path)

    # Initialize download statistics
    dl_stats = [0, 0]  # 0: success, 1: failed
    counter = 0

    # Determine the tile string (singular or plural)
    if len(search_result) == 1:
        tile_str = 'tile'
    else:
        tile_str = 'tiles'

    # Print the number of tiles to download
    print(f'\t- Downloading {len(search_result)} {tile_str}')
    print(f'\t\t- Requesting access token')
    token = get_access_token()

    # Iterate over the search results
    for key, value in search_result.items():
        # Create the save path and extract path
        save_path = os.path.join(base_path, 'tmp', f'{value["download_id"]}.zip')
        extract_path = os.path.join(base_path, f'R{value["relative_orbit"]:03}',
                                    f'{value["acquisition_date"].strftime("%Y%m%d")}')

        # Create the directories if they don't exist
        if not os.path.exists(os.path.dirname(save_path)):
            os.makedirs(os.path.dirname(save_path))
        if not os.path.exists(os.path.join(extract_path,
                                           f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')) or not os.path.exists(
            os.path.join(extract_path, f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')):
            # Print the scene name and progress
            print(f'\t\t- {value["scene_name"].replace(".SAFE", "")} ({counter + 1:0.0f}/{len(search_result):0.0f})')
            url = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products({value['download_id']})/$value"

            # Download the scene using wget
            checksum_check = False
            iteration = 0
            while not checksum_check and iteration <= 5:
                if iteration > 0:
                    print(f'\t\t\t- Checksums different. Requesting new token and redownloading (try {iteration}/5)')
                    token = get_access_token()
                iteration += 1
                os.system(f'wget -q --header "Authorization: Bearer {token}" \'{url}\' -O {save_path}')
                checksum_check = calculate_md5(save_path) == value["checksum"]

            if checksum_check:
                # Extract the required files
                tgt_filename_b4 = value["scene_name"].replace('.SAFE', '_B04_10m.jp2')
                tgt_filename_cloud = value["scene_name"].replace('.SAFE', '_MSK_CLDPRB_20m.jp2')

                print(f'\t\t\t- Extracting {tgt_filename_b4}')
                utils.extract_specific_files(save_path, extract_path,
                                             r'T\d{2}[A-Z]{3}_\d{4}\d{2}\d{2}T\d{2}\d{2}\d{2}_B04_10m\.jp2',
                                             f"{tgt_filename_b4}")
                print(f'\t\t\t- Extracting {tgt_filename_cloud}')
                utils.extract_specific_files(save_path, extract_path, r'MSK_CLDPRB_20m.jp2',
                                             f"{tgt_filename_cloud}")
                os.system(f'rm -rf {os.path.join(extract_path, value["scene_name"])}')
                os.system(f'rm -rf {save_path}')
                dl_stats[0] += 1
            else:
                print(
                    f'\t\t- {value["scene_name"].replace(".SAFE", "")} could not be downloaded successfully --> skipping ({counter + 1:0.0f}/{len(search_result):0.0f})')
                os.system(f'rm -rf {save_path}')
                dl_stats[1] += 1

        else:
            print(
                f'\t\t- {value["scene_name"].replace(".SAFE", "")} already present --> skipping ({counter + 1:0.0f}/{len(search_result):0.0f})')
            os.system(f'rm -rf {save_path}')

        counter += 1

    # Print the download statistics
    # if len(search_result)>0:
    #     print("\n".join([f'\t- Download statistics',f'\t\t- Successful: {dl_stats[0]}',f'\t\t- Failed: {dl_stats[1]}']))


def clean_search_result(search_result, base_path):
    """
    Cleans the search result by removing already downloaded tiles.

    Args:
        search_result (dict): The search result containing the scene information.
        base_path (str): The base path to check for already downloaded tiles.

    Returns:
        dict: The cleaned search result.
    """

    # Create a copy of the search result
    search_result_clean = search_result.copy()

    # Determine the tile string (singular or plural)
    if len(search_result_clean) == 1:
        tile_str = 'tile'
    else:
        tile_str = 'tiles'

    # Print the number of available tiles
    print(f'\t\t- {len(search_result_clean)} {tile_str} available')

    # Initialize a list to store the already downloaded tiles
    tiles_present = []

    # Iterate over the search result
    for key, value in search_result_clean.items():
        # Create the extract path
        extract_path = os.path.join(base_path, f'R{value["relative_orbit"]:03}',
                                    f'{value["acquisition_date"].strftime("%Y%m%d")}')

        # Check if the tile is already downloaded
        if os.path.exists(os.path.join(extract_path,
                                       f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')) and os.path.exists(
            os.path.join(extract_path, f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')):
            tiles_present.append(key)

    # Print the number of already downloaded tiles
    print(f'\t\t- {len(tiles_present)} {tile_str} already downloaded')

    # Remove the already downloaded tiles from the search result
    for key in tiles_present:
        search_result_clean.pop(key, 'None')

    # Return the cleaned search result
    return search_result_clean
