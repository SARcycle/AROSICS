import datetime
import glob
import os

import copernicus_api
import coreg_main
import utils

start_date = datetime.datetime(2024, 1, 1)
end_date = datetime.datetime.now()
mosaicing = False  # Set to False if calculation should be done on a tile-by-tile basis, True if it should be mosaiced beforehand
base_path = '/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/'
base_path_suffix = 'S2'

dates = [start_date + datetime.timedelta(days=i) for i in range((end_date - start_date).days + 1)]

for date in dates:
    str = f'{date.strftime("%Y-%m-%d %H:%M:%S")} - {(date + datetime.timedelta(days=1) - datetime.timedelta(seconds=1)).strftime("%Y-%m-%d %H:%M:%S")}'
    print("=" * len(str))
    print(str)
    print("-" * len(str))
    search_result = copernicus_api.image_search(date, date + datetime.timedelta(days=1))
    search_result_clean = copernicus_api.clean_search_result(search_result, os.path.join(base_path, base_path_suffix))
    if len(search_result_clean) > 0:
        copernicus_api.S2_scene_download(os.path.join(base_path, base_path_suffix), search_result=search_result_clean)

    unique_combinations = set()

    # Iterate over the dictionary
    for entry in search_result.values():
        # Add the combination of 'ac_date' and 'rel_robit' to the set
        unique_combinations.add((entry['acquisition_date'], entry['relative_orbit']))

    # Create a dictionary to store the counts
    counts = {}

    # Iterate over the dictionary again
    for entry in search_result.values():
        # Create a tuple for the current combination
        current_combination = (entry['acquisition_date'], entry['relative_orbit'])

        # If the combination is in the unique list, increment its count
        if current_combination in unique_combinations:
            if current_combination in counts:
                counts[current_combination] += 1
            else:
                counts[current_combination] = 1

    # Convert the set to a list
    unique_list = list(unique_combinations)
    for acquisitions in counts.items():
        relative_orbit = acquisitions[0][1]
        acquisition_date = acquisitions[0][0]
        tile_count = acquisitions[1]
        s2_bo4_tiles = glob.glob(os.path.join(os.path.join(base_path, base_path_suffix), f'R{relative_orbit:03}',
                                              f'{acquisition_date.strftime("%Y%m%d")}', '*_B04_10m.jp2'))
        s2_cloud_tiles = glob.glob(os.path.join(os.path.join(base_path, base_path_suffix), f'R{relative_orbit:03}',
                                                f'{acquisition_date.strftime("%Y%m%d")}', '*_MSK_CLDPRB_20m.jp2'))
        if len(s2_bo4_tiles) == tile_count and len(s2_cloud_tiles) == tile_count:  # If all files are present
            if mosaicing:
                b04_file_path_tmp, cld_file_path_tmp = coreg_main.S2_mosaic_scenes(
                    os.path.join(base_path, base_path_suffix), acquisition_date, relative_orbit)
                coreg_main.coregister_S2(b04_file_path_tmp, cld_file_path_tmp,
                                         os.path.join(os.path.join(base_path, base_path_suffix),
                                                      f'R{relative_orbit:03}',
                                                      f'{acquisition_date.strftime("%Y%m%d")}'), cloud_threshold=65,
                                         mosaicing=mosaicing)
            else:
                for idx in range(len(s2_bo4_tiles)):
                    if utils.get_epsg(s2_bo4_tiles[idx]) != 32632:
                        if not os.path.exists(s2_bo4_tiles[idx].replace(".jp2",
                                                                        "_utm32.tif")):  # Check if file is in UTM32 since result from GDALWARP with differing coordinate systems in input files is incorrect
                            os.system(
                                f'gdalwarp -tr 10 10 -r near -tap -q -ot UInt16 -s_srs EPSG:{utils.get_epsg(s2_bo4_tiles[idx])} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {s2_bo4_tiles[idx]} {s2_bo4_tiles[idx].replace(".jp2", "_utm32.tif")}')
                        b04_file_path_tmp = s2_bo4_tiles[idx].replace(".jp2", "_utm32.tif")
                    else:
                        b04_file_path_tmp = s2_bo4_tiles[idx]

                    if utils.get_epsg(s2_cloud_tiles[idx]) != 32632:
                        if not os.path.exists(s2_cloud_tiles[idx].replace("20m.jp2",
                                                                          "10m_utm32.tif")):  # Check if file is in UTM32 since result from GDALWARP with differing coordinate systems in input files is incorrect
                            os.system(
                                f'gdalwarp -tr 10 10 -r near -tap -q -ot Byte -s_srs EPSG:{utils.get_epsg(s2_cloud_tiles[idx])} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {s2_cloud_tiles[idx]} {s2_cloud_tiles[idx].replace("20m.jp2", "10m_utm32.tif")}')
                        cld_file_path_tmp = s2_cloud_tiles[idx].replace("20m.jp2", "10m_utm32.tif")
                    elif utils.get_epsg(s2_cloud_tiles[idx]) == 32632 and not os.path.exists(
                            s2_cloud_tiles[idx].replace("20m.jp2", "10m.tif")):
                        os.system(
                            f'gdalwarp -tr 10 10 -r near -q -ot Byte -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {s2_cloud_tiles[idx]} {s2_cloud_tiles[idx].replace("20m.jp2", "10m.tif")}')
                        cld_file_path_tmp = s2_cloud_tiles[idx].replace("20m.jp2", "10m.tif")
                    else:
                        cld_file_path_tmp = s2_cloud_tiles[idx].replace("20m.jp2", "10m.tif")

                    coreg_main.coregister_S2(b04_file_path_tmp, cld_file_path_tmp,
                                             os.path.join(os.path.join(base_path, base_path_suffix),
                                                          f'R{relative_orbit:03}',
                                                          f'{acquisition_date.strftime("%Y%m%d")}'), cloud_threshold=65,
                                             mosaicing=mosaicing)
