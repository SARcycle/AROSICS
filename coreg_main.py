# Import necessary libraries
import os

import numpy as np
from arosics import COREG_LOCAL
from osgeo import gdal

import utils


# Function to get the extent and dimensions of a raster

def S2_mosaic_scenes(base_path, acquisition_date, rel_orbit):
    import glob, os

    extract_path = os.path.join(base_path, f'R{rel_orbit:03}',
                                f'{acquisition_date.strftime("%Y%m%d")}')
    b04_path_out = os.path.join(extract_path, f'S2-L2A-mosaic_{acquisition_date.strftime("%Y-%m-%dT%H%M%S")}_B04-10m')
    cld_path_out = os.path.join(extract_path, f'S2-L2A-mosaic_{acquisition_date.strftime("%Y-%m-%dT%H%M%S")}_cloud-20m')

    if not os.path.exists(f'{b04_path_out}.tif'):  # If B04 mosaic missing
        # Mosaic of B04
        files_present_b04 = glob.glob(os.path.join(base_path, extract_path, '*_B04_10m.jp2'))
        files_present_b04_utm32 = []  # Create a list of files that are in UTM32
        for file in files_present_b04:
            if utils.get_epsg(
                    file) != 32632:  # Check if file is in UTM32 since result from GDALWARP with differing coordinate systems in input files is incorrect
                os.system(
                    f'gdalwarp -q -ot UInt16 -s_srs EPSG:{utils.get_epsg(file)} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {file} {file.replace(".jp2", "_utm32.tif")}')
                files_present_b04_utm32.append(file.replace('.jp2', '_utm32.tif'))
            else:
                files_present_b04_utm32.append(file)

        print(f'\t- Creating vrt ({os.path.basename(b04_path_out)}.vrt)')
        os.system(
            f'gdalbuildvrt -allow_projection_difference -q -srcnodata 0 {b04_path_out}.vrt {" ".join(files_present_b04_utm32)}')
        print(f'\t- Mosaiking ({os.path.basename(b04_path_out)}.tif)')
        os.system(
            f'gdalwarp -tr 10 10 -r near -tap -q -ot UInt16 -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {b04_path_out}.vrt {b04_path_out}.tif')
        os.system(f'rm {" ".join(glob.glob(os.path.join(base_path, extract_path, "*_utm32.tif")))}')
        # # Exchange in case COG is desired
        # print(f'\t- Mosaiking ({os.path.basename(b04_path_out)}.tif)')
        # os.system(f'gdalwarp -q -ot UInt16 -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=LZW -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES -co TILED=YES {b04_path_out}.vrt {b04_path_out}.tif')
        # print(f'\t- Converting to COG ({os.path.basename(b04_path_out)}.tif)')
        # os.system(f'gdal_translate -q -ot UInt16 -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=LZW -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS {b04_path_out}.vrt {b04_path_out}.tif')

    if not os.path.exists(f'{cld_path_out}.tif'):  # If cloud mosaic missing
        # Mosaic of Clouds
        files_present_cld = glob.glob(os.path.join(base_path, extract_path, '*_MSK_CLDPRB_20m.jp2'))
        files_present_cld_utm32 = []  # Create a list of files that are in UTM32
        for file in files_present_cld:
            if utils.get_epsg(
                    file) != 32632:  # Check if file is in UTM32 since result from GDALWARP with differing coordinate systems in input files is incorrect
                os.system(
                    f'gdalwarp -q -ot Byte -s_srs EPSG:{utils.get_epsg(file)} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {file} {file.replace(".jp2", "_utm32.tif")}')
                files_present_cld_utm32.append(file.replace('.jp2', '_utm32.tif'))
            else:
                files_present_cld_utm32.append(file)
        print(f'\t- Creating vrt ({os.path.basename(cld_path_out)}.vrt)')
        os.system(
            f'gdalbuildvrt -q -allow_projection_difference -q -srcnodata 0 {cld_path_out}.vrt {" ".join(files_present_cld_utm32)}')
        print(f'\t- Mosaiking ({os.path.basename(cld_path_out)}.tif)')
        os.system(
            f'gdalwarp -tr 10 10 -r near -tap -q -ot Byte -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS {cld_path_out}.vrt {cld_path_out}.tif')
        os.system(
            f'rm {" ".join(glob.glob(os.path.join(base_path, extract_path, "*_utm32.tif")))}')  # Deleting intermediate UTM32 files needed for mosaicing
    return f'{b04_path_out}.tif', f'{cld_path_out}.tif'


def get_extent_and_dimensions(filepath):
    # Open the raster file
    with gdal.Open(filepath, gdal.GA_ReadOnly) as ds:
        # Get the geotransform parameters
        geoTransform = ds.GetGeoTransform()
        # Calculate the extents
        minx = geoTransform[0]
        maxy = geoTransform[3]
        maxx = minx + geoTransform[1] * ds.RasterXSize
        miny = maxy + geoTransform[5] * ds.RasterYSize
        # Get the dimensions
        x_size, y_size = ds.RasterXSize, ds.RasterYSize
    # Return the extents and dimensions
    return (minx, maxx, miny, maxy, x_size, y_size)


# Function to equalize the extents of the reference and target images

def resample_raster(src_filename, match_filename, dst_filename):
    # Get information from the match file
    match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize
    match_ds = None

    # Create a temporary file to change the resolution
    temp_filename = "/tmp/temp.tif"
    res = match_geotrans[1]
    os.system(f'gdalwarp -tr {res} {res} -tap {src_filename} {temp_filename}')

    # Now cut the temporary file to match the extent of the match file
    os.system(
        f'gdalwarp -te {match_geotrans[0]} {match_geotrans[3] - high * res} {match_geotrans[0] + wide * res} {match_geotrans[3]} -dstnodata 0 {temp_filename} {dst_filename}')

    # Remove the temporary file
    os.system(f'rm {temp_filename}')


def equalize_extents(im_reference, im_target):
    # Get the extent and dimensions of the target image
    minx_target, maxx_target, miny_target, maxy_target, x_size_target, y_size_target = get_extent_and_dimensions(
        im_target)
    # Define the output file name
    outputFile = im_reference.replace('.tif', '_masked.tif')

    # If the output file already exists
    if os.path.exists(outputFile):
        # Get the extent and dimensions of the reference image
        minx_ref, maxx_ref, miny_ref, maxy_ref, x_size_ref, y_size_ref = get_extent_and_dimensions(im_reference)
        # If all extents and dimensions are equal, there's nothing to do
        if minx_target == minx_ref and maxx_target == maxx_ref and miny_target == miny_ref and maxy_target == maxy_ref and x_size_target == x_size_ref and y_size_target == y_size_ref:
            return outputFile

    # If the extents are not equal, crop the reference image to the target image size
    print('Cropping Reference Image to Target Image Size')
    kwargs = {"format": "GTiff", "outputBounds": [minx_target, miny_target, maxx_target, maxy_target]}
    with gdal.Warp(outputFile, im_reference, **kwargs) as ds:
        pass
    return outputFile


# Function to store X- and Y-Offsets as geotifs
def shifts_to_tif(CRL, outfolder=None):
    from osgeo import gdal
    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2']
    if outfolder is None:
        outfolder = os.path.dirname(CRL.path_out)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    if len(CRL.coreg_info['GCPList']) == 0:  # If there are no valid GCPs --> Clouds?
        with gdal.Open(CRL.im2shift.filePath) as ds:
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            x_shifts = np.full((rows, cols), 0)
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()
            driver = gdal.GetDriverByName('GTiff')
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dx.tif')),
                                  cols, rows, 1, gdal.GDT_Byte, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing X Shifts')
            band.WriteArray(x_shifts)

            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dy.tif')),
                                  cols, rows, 1, gdal.GDT_Byte, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing Y Shifts')  # Because both are zero-matrices
            band.WriteArray(x_shifts)
            del x_shifts
    else:
        print('Calculating X-Shifts')
        x_shifts = CRL.tiepoint_grid.to_interpolated_raster(metric='X_SHIFT_PX')
        with gdal.Open(CRL.im2shift.filePath) as ds:
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()
            driver = gdal.GetDriverByName('GTiff')
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dx.tif')),
                                  cols, rows, 1, gdal.GDT_Float32, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing X Shifts')
            band.WriteArray(x_shifts)
            del x_shifts
        print('Calculating Y-Shifts')
        y_shifts = CRL.tiepoint_grid.to_interpolated_raster(metric='Y_SHIFT_PX')
        with gdal.Open(CRL.im2shift.filePath) as ds:
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()
            driver = gdal.GetDriverByName('GTiff')
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dy.tif')),
                                  cols, rows, 1, gdal.GDT_Float32, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing Y Shifts')
            band.WriteArray(y_shifts)
            del y_shifts


def coregister_S2(file_path_in, cld_mask_path_in, folder_path_out, cloud_threshold=65, mosaicing=False):
    max_points = 5000
    im_target = file_path_in
    if mosaicing:
        cld_mask_path_out = cld_mask_path_in.replace(".jp2", "_bin.tif")
    else:
        cld_mask_path_out = cld_mask_path_in.replace(".tif", "_bin.tif")

    x_size, y_size = get_extent_and_dimensions(im_target)[-2:]

    grid_res = round((((x_size * y_size) / max_points) ** .5) / 5) * 5
    if mosaicing:
        out_name = os.path.basename(file_path_in).replace('B04', 'registration_swiss')
    else:
        out_name = os.path.basename(file_path_in).replace('B04', 'registration_swiss').replace('jp2', 'tif')
    out_folder = os.path.join(folder_path_out, f'TiePoint_GridRes_{grid_res}x{grid_res}px')
    if not os.path.exists(os.path.join(out_folder, out_folder)) or not os.path.exists(
            os.path.join(out_folder, out_name.replace(".tif", ".shp"))) or not os.path.exists(
        os.path.join(out_folder, out_name.replace(".tif", "_dx.tif"))) or not os.path.exists(
        os.path.join(out_folder, out_name.replace(".tif", "_dy.tif"))):

        # print(f'\t- Converting Cloud Probability to Cloud Mask (≥{cloud_threshold}%): ', end="")
        # os.system(f'gdal_calc.py -A {cld_mask_path_in} --overwrite --outfile={cld_mask_path_out} --calc="logical_and(A>={cloud_threshold}, A<=100)" --type=Byte --NoDataValue=0 --co COMPRESS=LZW')
        im_reference = '/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/S2_GRI.tif'
        # Calculate the grid resolution based on a maximum point count and the raster dimensions

        # Print the title
        title_str = f'{os.path.basename(im_target)} | {grid_res}x{grid_res}px'
        print('=' * len(title_str))
        print(title_str)
        print('-' * len(title_str))

        # Equalize the extents of the reference and target images
        im_reference_masked = equalize_extents(im_reference, im_target)

        # Define the arguments for the coregistration
        kwargs = {
            'grid_res': grid_res,
            'path_out': f'{out_name}',
            'projectDir': out_folder,
            'q': False,
            'nodata': (0, 0),
            # 'mask_baddata_tgt': cld_mask_path_out,
            'CPUs': max(os.cpu_count() - 1, 1),
            # Automatically use max number of cores, but leave one core for system (if more than 1 core CPU)
            'progress': False,
        }

        # Perform the coregistration
        CRL = COREG_LOCAL(im_reference_masked, im_target, **kwargs)
        if len(CRL.coreg_info['GCPList']) != 0:  # If there are GCPs, so not totally cloud covered
            CRL.correct_shifts()
            CRL.tiepoint_grid.to_PointShapefile(
                path_out=os.path.join(folder_path_out, f'TiePoint_GridRes_{grid_res}x{grid_res}px',
                                      out_name.replace(".tif", ".shp")))

        shifts_to_tif(CRL)

        # # AWS Cost Estimation
        # # https://calculator.aws/#/estimate?id=0a2b8d3c19078290f471fbba9938fa6e2235e875
        print('=' * len(title_str))

    # start = time.time()
    # # Define the location and base path
    # location = 'CH'
    # base_path = '/mnt/c/Users/Localadmin/Documents/SATROMO/AROSICS_Coregistration/Test/'
    # max_points = 5000
    # min_rel = 60
    # im_reference = '/mnt/c/Users/Localadmin/Documents/SATROMO/AROSICS_Coregistration/Test/CH/input/S2_GRI.tif'
    #
    # # Iterate over all files with the corresponding naming scheme
    # for im_target in glob.glob(base_path + location + '/input/20200723*B4_10m_mosaic.tif'):
    #     # Get the dimensions of the target image
    #     x_size, y_size = get_extent_and_dimensions(im_target)[-2:]
    #     # Calculate the grid resolution based on a maximum point count and the raster dimensions
    #     grid_res = round((((
    #                                    x_size * y_size) / max_points) ** .5) / 5) * 5
    #
    #     # Print the title
    #     title_str = f'{os.path.basename(im_target)} | {grid_res}x{grid_res}px | {min_rel}Pct_minReliability'
    #     print('=' * len(title_str))
    #     print(title_str)
    #     print('-' * len(title_str))
    #
    #     # Equalize the extents of the reference and target images
    #     im_reference_masked = equalize_extents(im_reference, im_target)
    #
    #     # Define the arguments for the coregistration
    #     kwargs = {
    #         'grid_res': grid_res,
    #         'path_out': 'auto',
    #         'projectDir': base_path + location + '/output/TiePoint_GridRes_' + str(grid_res) + 'x' + str(
    #             grid_res) + 'px_' + str(min_rel) + 'Pct_minReliability',
    #         'q': False,
    #         'nodata': (0, 0),
    #     }
    #
    #     # Perform the coregistration
    #     CRL = COREG_LOCAL(im_reference_masked, im_target, **kwargs)
    #     print(CRL)
    #     CRL.correct_shifts()
    #
    #     CRL.tiepoint_grid.to_PointShapefile(
    #         path_out=base_path + location + '/output/TiePoint_GridRes_' + str(grid_res) + 'x' + str(grid_res) + 'px_' + str(
    #             min_rel) + 'Pct_minReliability/' + os.path.basename(im_target).replace('tif', 'shp'))
    #
    #     shifts_to_tif(CRL)
    #
    #     # # AWS Cost Estimation
    #     # # https://calculator.aws/#/estimate?id=0a2b8d3c19078290f471fbba9938fa6e2235e875
    #     print('=' * len(title_str))
    #     print(f'Completed in {time.time()-start:0.3f} sec')
    else:
        print('\t- Image already coregistered --> Skipping')


def mosaic_tiles_old(input_folder, operation, bands):
    import os
    import datetime
    import rasterio
    import rioxarray

    import warnings
    warnings.filterwarnings('ignore')

    # Get a list of all tile files
    tile_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(('_dx.tif', '_dy.tif')):
                tile_files.append(os.path.join(root, file))

    # Separate dx and dy files
    dx_files = [file for file in tile_files if file.endswith('_dx.tif')]
    dy_files = [file for file in tile_files if file.endswith('_dy.tif')]

    # Create a dictionary to store the mosaic data

    # Process dx files
    if bands in ['dx', 'both']:
        dx_datasets = []
        for file in dx_files:
            dx_datasets.append(rioxarray.open_rasterio(file, mask_and_scale=True))
        if dx_datasets:
            dx_sum = rioxarray.merge.merge_arrays(dataarrays=dx_datasets, method='sum', nodata=0)
            dx_count = rioxarray.merge.merge_arrays(dataarrays=dx_datasets, method='count', nodata=0)
            dx_avg = dx_sum / dx_count
            del dx_sum, dx_count

    # Process dy files
    if bands in ['dy', 'both']:
        dy_datasets = []
        for file in dy_files:
            dy_datasets.append(rioxarray.open_rasterio(file, mask_and_scale=True))
        if dy_datasets:
            dy_sum = rioxarray.merge.merge_arrays(dataarrays=dy_datasets, method='sum', nodata=0)
            dy_count = rioxarray.merge.merge_arrays(dataarrays=dy_datasets, method='count', nodata=0)
            dy_avg = dy_sum / dy_count
            del dy_sum, dy_count

    # Create the output file name
    output_folder = input_folder
    if len(dx_files) > 0:
        dx_basename = os.path.basename(dx_files[1])
        datestr_dt = datetime.datetime.strptime(dx_basename.split('_')[2],
                                                '%Y%m%dT%H%M%S')  # Extracting the time value of the input file
        datestr = datestr_dt.strftime('%Y-%m-%dt%H%M%S')
    else:
        dy_basename = os.path.basename(dy_files[1])
        datestr_dt = datetime.datetime.strptime(dy_basename.split('_')[2],
                                                '%Y%m%dT%H%M%S')  # Extracting the time value of the input file
        datestr = datestr_dt.strftime('%Y-%m-%dt%H%M%S')

    if bands == 'both':
        output_filename = f'S2-L2A-mosaic_{datestr}_registration_swiss-10m.tif'
    elif bands == 'dx':
        output_filename = f'S2-L2A-mosaic_{datestr}_registration_swiss-10m_dx.tif'
    elif bands == 'dy':
        output_filename = f'S2-L2A-mosaic_{datestr}_registration_swiss-10m_dy.tif'

    # Write the mosaic to a new file
    if bands == 'both':
        dy_avg.rio.to_raster("rx_merged_GLC_FCS30_ave.tiff")
    else:
        with rasterio.open(os.path.join(output_folder, output_filename), 'w', driver='GTiff',
                           height=mosaic_data[bands].shape[0], width=mosaic_data[bands].shape[1], count=1,
                           dtype=rasterio.float32, crs='EPSG:32632',
                           transform=rasterio.transform.from_bounds(0, 0, 1000, 1000, 1000, 1000)) as dst:
            dst.write(mosaic_data[bands], 1)

    # IDea

    import gdal
    import numpy as np

    # Merge the tiles using gdal_merge
    gdal_merge_str = "gdal_merge.py -o merged.tif -of GTiff *.tif"
    subprocess.call(gdal_merge_str, shell=True)

    # Open the merged raster
    ds = gdal.Open("merged.tif")
    band = ds.GetRasterBand(1)
    width, height = ds.RasterXSize, ds.RasterYSize

    # Create an ndarray to store the summed values and counters
    sum_array = np.zeros((2, height, width), dtype=np.float32)

    # Iterate over the tiles
    for tile in ["tile1.tif", "tile2.tif", ...]:  # list of tile files
        tile_ds = gdal.Open(tile)
        tile_band = tile_ds.GetRasterBand(1)
        tile_data = tile_band.ReadAsArray()

        # Check if the tile data type is not byte
        if tile_band.DataType != gdal.GDT_Byte:
            # Add the tile values to the sum array
            sum_array[0] += tile_data
            # Increment the counter array
            sum_array[1] += (tile_data != 0).astype(np.int8)

    # Calculate the average values
    avg_array = sum_array[0] / sum_array[1]

    # Save the average array to a new GeoTIFF file
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.CreateCopy("average.tif", ds, 1, [avg_array])
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())


def mosaic_tiles(input_folder, operation, bands):
    # Creating reference mosaic from all input files
    import datetime

    # Get a list of all tile files
    tile_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(('_dx.tif', '_dy.tif')):
                tile_files.append(os.path.join(root, file))

    # Separate dx and dy files
    dx_files = [file for file in tile_files if file.endswith('_dx.tif')]
    dy_files = [file for file in tile_files if file.endswith('_dy.tif')]

    # Create the output file name
    output_folder = input_folder
    if len(dx_files) > 0:
        dx_basename = os.path.basename(dx_files[1])
        datestr_dt = datetime.datetime.strptime(dx_basename.split('_')[2],
                                                '%Y%m%dT%H%M%S')  # Extracting the time value of the input file
        datestr = datestr_dt.strftime('%Y-%m-%dt%H%M%S')
    else:
        dy_basename = os.path.basename(dy_files[1])
        datestr_dt = datetime.datetime.strptime(dy_basename.split('_')[2],
                                                '%Y%m%dT%H%M%S')  # Extracting the time value of the input file
        datestr = datestr_dt.strftime('%Y-%m-%dt%H%M%S')

    if bands == 'both':
        output_filename = f'S2-L2A-mosaic_{datestr}_registration_swiss-10m.tif'
    elif bands == 'dx':
        output_filename = f'S2-L2A-mosaic_{datestr}_registration_swiss-10m_dx.tif'
    elif bands == 'dy':
        output_filename = f'S2-L2A-mosaic_{datestr}_registration_swiss-10m_dy.tif'

    # Create a reference image with the correct dimensions, data type and coordinate axis
    os.system(
        f'gdal_merge.py -init 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -o {os.path.join(output_folder, output_filename)} -ot Int32 -createonly {" ".join(dx_files)}')

    combine_geotiffs(output_folder=output_folder, output_filename=output_filename, input_geotiff_paths=dx_files,
                     operation='average')

    pass


def combine_geotiffs(output_folder, output_filename, input_geotiff_paths, operation):
    import os
    import numpy as np
    from osgeo import gdal

    # Open the output GeoTIFF
    output_ds = gdal.Open(os.path.join(output_folder, output_filename), gdal.GA_Update)
    output_band = output_ds.GetRasterBand(1)

    # Create a numpy array based on the metadata
    output_data_type = get_gdal_raster_data_types(output_band.DataType)
    output_array = np.zeros((output_ds.RasterYSize, output_ds.RasterXSize), dtype=output_data_type)

    # Initialize a count array to keep track of the number of valid values
    count_array = np.zeros((output_ds.RasterYSize, output_ds.RasterXSize), dtype=np.int16)

    # Loop through the input GeoTIFFs
    for input_geotiff_path in input_geotiff_paths:
        # Open the input GeoTIFF
        input_ds = gdal.Open(input_geotiff_path)
        input_band = input_ds.GetRasterBand(1)

        input_data_type = get_gdal_raster_data_types(input_band.DataType)

        # Check if the input data type is Byte and the value is 0 (no data)
        if input_data_type == np.byte and np.all(np.array(input_band.ReadAsArray()) == 0):
            continue

        # Get the geotransform from the input GeoTIFF
        input_gt = input_ds.GetGeoTransform()

        # Calculate the offset to align the input data with the output data
        offset_x = int((input_gt[0] - output_ds.GetGeoTransform()[0]) / output_ds.GetGeoTransform()[1])
        offset_y = int((input_gt[3] - output_ds.GetGeoTransform()[3]) / output_ds.GetGeoTransform()[5])

        # Check if the input data is within the bounds of the output data
        if offset_x < 0 or offset_y < 0 or offset_x + input_band.XSize > output_ds.RasterXSize or offset_y + input_band.YSize > output_ds.RasterYSize:
            raise ValueError("Input data is outside the bounds of the output data")

        # Calculate the chunk size
        chunk_size = 1024
        for i in range(0, input_band.YSize, chunk_size):
            for j in range(0, input_band.XSize, chunk_size):
                # Calculate the chunk dimensions
                chunk_x = min(chunk_size, input_band.XSize - j)
                chunk_y = min(chunk_size, input_band.YSize - i)

                # Read the input chunk
                chunk = input_band.ReadAsArray(j, i, chunk_x, chunk_y)
                # Conversion from pixels to shift in cm based on pixel size
                chunk = chunk * input_gt[1] * 100

                if operation == 'average':
                    output_array[offset_y + i:offset_y + i + chunk_y,
                    offset_x + j:offset_x + j + chunk_x] += output_data_type(chunk)
                    count_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x] += 1
                elif operation == 'median':
                    output_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x] = np.median(
                        np.dstack((output_array[offset_y + i:offset_y + i + chunk_y,
                                   offset_x + j:offset_x + j + chunk_x], chunk)), axis=2)[0]
                elif operation == 'sum':
                    output_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x] += chunk
                elif operation == 'max':
                    output_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x] = np.maximum(
                        output_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x], chunk)
                elif operation == 'min':
                    output_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x] = np.minimum(
                        output_array[offset_y + i:offset_y + i + chunk_y, offset_x + j:offset_x + j + chunk_x], chunk)
                else:
                    raise ValueError("Invalid operation. Must be one of 'average', 'median', 'sum', 'max', or 'min'")

    # If the operation is 'average', divide by the count array
    if operation == 'average':
        output_array[count_array > 0] = output_array[count_array > 0] / count_array[count_array > 0]

    # Set cells without a valid value to 0
    output_array[count_array == 0] = 0

    # Write the combined array back to the GeoTIFF
    output_band.WriteArray(output_array)

    # Close the datasets
    output_ds = None
    for input_ds in [gdal.Open(path) for path in input_geotiff_paths]:
        input_ds = None


def get_gdal_raster_data_types(gdal_data_type):
    # Finding the output data type in accordance with http://ikcest-drr.osgeo.cn/tutorial/k8023
    if gdal_data_type == 1:
        return np.byte
    elif gdal_data_type == 2:
        return np.uint16
    elif gdal_data_type == 3:
        return np.int16
    elif gdal_data_type == 4:
        return np.uint32
    elif gdal_data_type == 5:
        return np.int32
    elif gdal_data_type == 6:
        return np.float32
    elif gdal_data_type == 7:
        return np.float64
    else:
        return None


def calculate_idw(shapefile_folder, reference_geotiff, output_raster, power=2, max_distance=None, max_points=None,
                  min_reliability=60):
    import os
    import numpy as np
    from osgeo import ogr
    import rasterio, time

    with rasterio.open(reference_geotiff) as src:
        transform = src.transform
        rows, cols = src.height, src.width
        x_coords_px, y_coords_px = rasterio.transform.xy(transform, range(cols),
                                                         range(rows))  # Referenced to pixel center

    x_shift_array = np.zeros((rows, cols), dtype=np.int32)
    y_shift_array = np.zeros((rows, cols), dtype=np.int32)

    # Merge shapefiles from subfolders
    merged_features = []
    for root, dirs, files in os.walk(shapefile_folder):
        for filename in files:
            if filename.endswith(".shp"):
                shapefile_path = os.path.join(root, filename)
                shape_ds = ogr.Open(shapefile_path)
                layer = shape_ds.GetLayer()
                merged_features.extend(layer)

    filtered_features = [feature for feature in merged_features if
                         feature.GetField("OUTLIER") != 1 and feature.GetField("RELIABILIT") >= min_reliability]
    del merged_features

    x_coords_pts = np.array([feature.GetGeometryRef().GetX() for feature in filtered_features])
    y_coords_pts = np.array([feature.GetGeometryRef().GetY() for feature in filtered_features])
    x_shift_pts = np.array([feature.GetField('X_SHIFT_M') for feature in filtered_features])
    y_shift_pts = np.array([feature.GetField('Y_SHIFT_M') for feature in filtered_features])

    start = time.time()
    for x_count, x in enumerate(x_coords_px):
        if x_count > 0:
            print(
                f'{x_count / len(x_coords_px)} in {time.time() - start} sec (--> {((time.time() - start) / x_count) * len(x_coords_px)} remaining)')

        for y_count, y in enumerate(y_coords_px):
            x_mask = np.abs(x - x_coords_pts) <= max_distance
            y_mask = np.abs(y - y_coords_pts) <= max_distance
            indices = np.where(x_mask & y_mask)[0]
            if len(indices) > 0:
                x_coords_pts_tmp = x_coords_pts[indices]
                y_coords_pts_tmp = y_coords_pts[indices]
                x_shift_pts_tmp = x_shift_pts[indices]
                y_shift_pts_tmp = y_shift_pts[indices]

                # Calculate IDW value and select closest number of points (n=max_points)
                distances_tmp = np.sqrt((x - x_coords_pts_tmp) ** 2 + (y - y_coords_pts_tmp) ** 2)
                if max_points is not None:
                    sorted_indices = np.argsort(distances_tmp)
                    point_slct = sorted_indices[:max_points]
                    distances_tmp = distances_tmp[point_slct]
                    x_shift_pts_tmp = x_shift_pts[point_slct]
                    y_shift_pts_tmp = y_shift_pts[point_slct]
                weights = 1 / (distances_tmp ** power)
                total_weights = np.sum(weights)
                x_shift_array[y_count, x_count] = np.int32((np.sum(x_shift_pts_tmp * weights) / total_weights) * 100)
                y_shift_array[y_count, x_count] = np.int32((np.sum(y_shift_pts_tmp * weights) / total_weights) * 100)
            else:
                x_shift_array[y_count, x_count] = 0
                y_shift_array[y_count, x_count] = 0

    # Write IDW value to output raster
    # output_ds.GetRasterBand(1).WriteArray(idw_value, xoff=int(x), yoff=int(y))

    output_ds = None
    print(f"IDW interpolation saved to {output_raster}")
