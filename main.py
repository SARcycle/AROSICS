import datetime
import glob
import os
import time

import copernicus_api
import coreg_main
import utils

# Define start and end dates for the data processing
start_date = datetime.datetime(2021, 1, 1)
end_date = datetime.datetime(2021, 7, 31)
# end_date = datetime.datetime.now()

# Set mosaicing flag to False for tile-by-tile processing or True for mosaicing
mosaicing = True

# Define base path and suffix for the data
base_path = '/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/'
base_path_suffix = 'S2'

# Create a list of dates to process
dates = [start_date + datetime.timedelta(days=i) for i in range((end_date - start_date).days + 1)]

# Iterate over the dates
for date in dates:
    try:
        start = time.time()
        # Print a header for the current date
        str = f'{date.strftime("%Y-%m-%d %H:%M:%S")} - {(date + datetime.timedelta(days=1) - datetime.timedelta(seconds=1)).strftime("%Y-%m-%d %H:%M:%S")}'
        print("=" * len(str))
        print(str)
        print("-" * len(str))

        # Search for Copernicus API images for the current date
        search_result = copernicus_api.image_search(date, date + datetime.timedelta(days=1))
        search_result_clean = copernicus_api.clean_search_result(search_result,
                                                                 os.path.join(base_path, base_path_suffix))

        # Download the images if the search result is not empty
        if len(search_result_clean) > 0:
            copernicus_api.S2_scene_download(os.path.join(base_path, base_path_suffix),
                                             search_result=search_result_clean)

        # Create a set to store unique combinations of acquisition date and relative orbit
        unique_combinations = set()

        # Iterate over the search result dictionary
        for entry in search_result.values():
            # Add the combination of 'acquisition_date' and 'relative_orbit' to the set
            unique_combinations.add((entry['acquisition_date'], entry['relative_orbit']))

        # Create a dictionary to store the counts of each combination
        counts = {}

        # Iterate over the search result dictionary again
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

        # Iterate over the unique combinations and their counts
        for acquisitions in counts.items():
            relative_orbit = acquisitions[0][1]
            acquisition_date = acquisitions[0][0]
            tile_count = acquisitions[1]

            # Find the S2 BO4 and cloud tiles for the current combination
            s2_b04_tiles = glob.glob(os.path.join(os.path.join(base_path, base_path_suffix), f'R{relative_orbit:03}',
                                                  f'{acquisition_date.strftime("%Y%m%d")}', '*_B04_10m.jp2'))
            s2_cloud_tiles = glob.glob(os.path.join(os.path.join(base_path, base_path_suffix), f'R{relative_orbit:03}',
                                                    f'{acquisition_date.strftime("%Y%m%d")}', '*_MSK_CLDPRB_20m.jp2'))

            # Check if all files are present
            if len(s2_b04_tiles) == tile_count and len(s2_cloud_tiles) == tile_count:
                if mosaicing:
                    # Mosaic the scenes
                    b04_file_path_tmp, cld_file_path_tmp = coreg_main.S2_mosaic_scenes(
                        os.path.join(base_path, base_path_suffix), acquisition_date, relative_orbit)
                    coreg_main.coregister_S2(b04_file_path_tmp, cld_file_path_tmp,
                                             os.path.join(os.path.join(base_path, base_path_suffix),
                                                          f'R{relative_orbit:03}',
                                                          f'{acquisition_date.strftime("%Y%m%d")}'), cloud_threshold=65,
                                             mosaicing=mosaicing)
                else:
                    # Process each tile individually
                    for idx in range(len(s2_b04_tiles)):
                        # Check if the file needs to be reprojected to UTM32
                        if utils.get_epsg(s2_b04_tiles[idx]) != 32632:
                            if not os.path.exists(s2_b04_tiles[idx].replace(".jp2",
                                                                            "_utm32.tif")):
                                pixel_spacing_x, pixel_spacing_y = utils.get_pixel_spacing(s2_b04_tiles[idx])
                                os.system(
                                    f'gdalwarp -tr {pixel_spacing_x} {pixel_spacing_y} -overwrite  -r near -tap -q -ot UInt16 -s_srs EPSG:{utils.get_epsg(s2_b04_tiles[idx])} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {s2_b04_tiles[idx]} {s2_b04_tiles[idx].replace(".jp2", "_utm32.tif")}')
                            b04_file_path_tmp = s2_b04_tiles[idx].replace(".jp2", "_utm32.tif")
                        else:
                            b04_file_path_tmp = s2_b04_tiles[idx]

                        # Check if the cloud file needs to be reprojected to UTM32
                        if utils.get_epsg(s2_cloud_tiles[idx]) != 32632:
                            if not os.path.exists(s2_cloud_tiles[idx].replace("20m.jp2",
                                                                              "10m_utm32.tif")):
                                pixel_spacing_x, pixel_spacing_y = utils.get_pixel_spacing(s2_b04_tiles[idx])
                                os.system(
                                    f'gdalwarp -tr {pixel_spacing_x} {pixel_spacing_y} -r near -tap -q -ot Byte -overwrite -s_srs EPSG:{utils.get_epsg(s2_cloud_tiles[idx])} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {s2_cloud_tiles[idx]} {s2_cloud_tiles[idx].replace("20m.jp2", "10m_utm32.tif")}')
                            cld_file_path_tmp = s2_cloud_tiles[idx].replace("20m.jp2", "10m_utm32.tif")
                        elif utils.get_epsg(s2_cloud_tiles[idx]) == 32632 and not os.path.exists(
                                s2_cloud_tiles[idx].replace("20m.jp2", "10m.tif")):
                            pixel_spacing_x, pixel_spacing_y = utils.get_pixel_spacing(s2_b04_tiles[idx])
                            os.system(
                                f'gdalwarp -tap -tr {pixel_spacing_x} {pixel_spacing_y} -r near -q -ot Byte -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {s2_cloud_tiles[idx]} {s2_cloud_tiles[idx].replace("20m.jp2", "10m.tif")}')
                            cld_file_path_tmp = s2_cloud_tiles[idx].replace("20m.jp2", f"10m.tif")
                        else:
                            cld_file_path_tmp = s2_cloud_tiles[idx].replace("20m.jp2", f"10m.tif")

                        # Coregister the S2 BO4 and cloud tiles
                        coreg_main.coregister_S2(b04_file_path_tmp, cld_file_path_tmp,
                                                 os.path.join(os.path.join(base_path, base_path_suffix),
                                                              f'R{relative_orbit:03}',
                                                              f'{acquisition_date.strftime("%Y%m%d")}'),
                                                 cloud_threshold=65,
                                                 mosaicing=mosaicing)
        print(
            f'Total processing time for {datetime.datetime.strftime(date, "%Y/%m/%d")}: {datetime.timedelta(seconds=round(time.time() - start))}')
    except:
        print(f'!!!!!!!!!!!! Error on {datetime.datetime.strftime(date, "%Y/%m/%d")} !!!!!!!!!!!!')
