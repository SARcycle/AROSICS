import os.path
import sys

import utils


def image_search(start_date, end_date):
    import geopandas as gpd
    import requests
    from datetime import datetime
    print(f'\t- Searching for images')
    if isinstance(start_date, datetime):  # Check if date is a datetime object
        start_date = start_date.strftime('%Y-%m-%d')
    if isinstance(end_date, datetime):  # Check if date is a datetime object
        end_date = end_date.strftime('%Y-%m-%d')
    proc_lvl = 'L2A'
    data_collection = "SENTINEL-2"
    gdf = gpd.read_file(
        '/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/swissBOUNDARIES3D_1_5_LV95_LN02_buffered_5000m_simplified_DP_5000m.gpkg',
        layer='bbox')
    aoi = gdf['geometry'][0]
    request_url = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Collection/Name eq '{data_collection}' and OData.CSC.Intersects(area=geography'SRID=4326;{aoi}') and ContentDate/Start gt {start_date}T00:00:00.000Z and ContentDate/Start lt {end_date}T00:00:00.000Z&$top=1000"
    json = requests.get(request_url).json()
    ids = [d['Id'] for d in json['value'] if 'Id' in d and proc_lvl in d.get('Name', '')]
    acquisition_dates = [datetime.strptime(d['ContentDate']['Start'], '%Y-%m-%dT%H:%M:%S.%fZ') for d in json['value'] if
                         'ContentDate' in d and proc_lvl in d.get('Name', '')]
    scene_names = [d['Name'] for d in json['value'] if 'Name' in d and proc_lvl in d.get('Name', '')]
    checksums = [d['Checksum'][0]['Value'] for d in json['value'] if
                 'Checksum' in d and proc_lvl in d.get('Name', '')]  # MD5 checksum
    rel_orbits = []
    for scene_name in scene_names:
        parts = scene_name.split('_')
        if len(parts) >= 5:
            rel_orbits.append(int(parts[4].lstrip("R")))

    # Creating a dictionary
    search_result = {}
    for i in range(len(ids)):
        search_result[f'tile_{i:04}'] = {
            'download_id': ids[i],
            'acquisition_date': acquisition_dates[i],
            'relative_orbit': rel_orbits[i],
            'scene_name': scene_names[i],
            'checksum': checksums[i],
        }

    return search_result


def get_access_token():
    import utils, requests
    username, password = utils.get_credentials('copernicus')
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
            f"Access token creation failed. Reponse from the server was: {r.json()}"
        )
    return r.json()["access_token"]


def calculate_md5(file_path):
    import hashlib
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def S2_scene_download(base_path, start_date=None, end_date=None, search_result=None):
    import requests
    from datetime import datetime
    from requests.packages.urllib3.exceptions import InsecureRequestWarning

    if search_result is not None:
        pass
    elif start_date is None:
        print('No start date provided. Exiting.')
        sys.exit()
    elif end_date is None:
        end_date = datetime.datetime.now()
        print('No end date provided. Defaulting to current date.')
    else:
        requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
        if isinstance(start_date, datetime):  # Check if date is a datetime object
            start_date = start_date.strftime('%Y-%m-%d')
        if isinstance(end_date, datetime):  # Check if date is a datetime object
            end_date = end_date.strftime('%Y-%m-%d')

        search_result = image_search(start_date=start_date, end_date=end_date)
        search_result = clean_search_result(search_result, base_path)

    dl_stats = [0, 0]  # Statistics for the downloads, 1: success, 2: failed
    counter = 0

    if len(search_result) == 1:
        tile_str = 'tile'
    else:
        tile_str = 'tiles'
    print(f'\t- Downloading {len(search_result)} {tile_str}')
    print(f'\t\t- Requesting access token')
    token = get_access_token()

    for key, value in search_result.items():
        save_path = os.path.join(base_path, 'tmp', f'{value["download_id"]}.zip')
        extract_path = os.path.join(base_path, f'R{value["relative_orbit"]:03}',
                                    f'{value["acquisition_date"].strftime("%Y%m%d")}')

        if not os.path.exists(os.path.dirname(save_path)):
            os.makedirs(os.path.dirname(save_path))
        if not os.path.exists(os.path.join(extract_path,
                                           f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')) or not os.path.exists(
            os.path.join(extract_path, f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')):
            # session = requests.Session()
            # session.headers.update({'Authorization': f'Bearer {token}'})
            print(f'\t\t- {value["scene_name"].replace(".SAFE", "")} ({counter + 1:0.0f}/{len(search_result):0.0f})')
            url = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products({value['download_id']})/$value"
            # response = session.get(url, allow_redirects=False)
            # while response.status_code in (301, 302, 303, 307):
            #     url = response.headers['Location']
            #     print("\t- Redirected to:", url)
            #     response = session.get(url, allow_redirects=False, timeout=5, verify=False)
            #
            # # Create a progress bar for the download
            # total_size = int(response.headers.get("content-length", 0))
            # block_size = 1024
            #
            # # Use tqdm to display the progress
            # with tqdm(total=total_size, unit="B", unit_scale=True) as progress_bar:
            #     with open(save_path, "wb") as p:
            #         for data in response.iter_content(block_size):
            #             progress_bar.update(len(data))  # Update the progress bar
            #             p.write(data)
            #
            # # Verify that the entire file was downloaded
            # if total_size != 0 and progress_bar.n != total_size:
            #     raise RuntimeError("Could not download the file completely")
            checksum_check = False
            iteration = 0
            while not checksum_check and iteration <= 5:
                if iteration > 0:
                    print(f'\t\t\t- Checksums different. Requesting new token and redownloading (try {iteration}/5)')
                    token = get_access_token()
                iteration += 1
                os.system(f'wget -q --header "Authorization: Bearer {token}" \'{url}\' -O {save_path}')
                # response = requests.get(url, headers={"Authorization": f"Bearer {token}"})
                # with open(save_path, 'wb') as file:
                #     file.write(response.content)
                checksum_check = calculate_md5(save_path) == value["checksum"]

            if checksum_check:
                # tgt_filename_xml = value["scene_name"].replace('.SAFE', '_L2A_quality.xml')
                tgt_filename_b4 = value["scene_name"].replace('.SAFE', '_B04_10m.jp2')
                tgt_filename_cloud = value["scene_name"].replace('.SAFE', '_MSK_CLDPRB_20m.jp2')

                # print(f'\t\t- Extracting {tgt_filename_xml}')
                # utils.extract_specific_files(save_path, extract_path, r'L2A_QUALITY.xml',
                #                             f"{tgt_filename_xml}")

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
                # coreg_main.coregister_S2(os.path.join(extract_path, f"{tgt_filename_b4}"), os.path.join(extract_path, f"{tgt_filename_cloud}"), extract_path, 65)
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


#   if len(search_result)>0:
#       print("\n".join([f'\t- Download statistics',f'\t\t- Successful: {dl_stats[0]}',f'\t\t- Failed: {dl_stats[1]}']))


def clean_search_result(search_result, base_path):
    search_result_clean = search_result.copy()

    # Deleting already downloaded tiles from download list
    if len(search_result_clean) == 1:
        tile_str = 'tile'
    else:
        tile_str = 'tiles'

    print(f'\t\t- {len(search_result_clean)} {tile_str} available')
    tiles_present = []
    for key, value in search_result_clean.items():
        extract_path = os.path.join(base_path, f'R{value["relative_orbit"]:03}',
                                    f'{value["acquisition_date"].strftime("%Y%m%d")}')
        if os.path.exists(os.path.join(extract_path,
                                       f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')) and os.path.exists(
            os.path.join(extract_path, f'{value["scene_name"].replace(".SAFE", "_B04_10m.jp2")}')):
            tiles_present.append(key)
    print(f'\t\t- {len(tiles_present)} {tile_str} already downloaded')
    for key in tiles_present:
        search_result_clean.pop(key, 'None')

    return search_result_clean
# S2_scene_download('2024-04-01', '2024-06-05')
# S2_mosaic_scenes('/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2', acquisition_date, rel_orbit)
