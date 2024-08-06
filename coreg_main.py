# Import necessary libraries
import glob
import os

import numpy as np
from arosics import COREG_LOCAL
from osgeo import gdal, gdalconst
from scipy.interpolate import interp2d

import utils


# Function to get the extent and dimensions of a raster
def S2_mosaic_scenes(base_path, acquisition_date, rel_orbit):
    """
    This function mosaics Sentinel-2 scenes for a given base path, acquisition date, and relative orbit.

    Parameters:
    base_path (str): The base path where the scenes are located.
    acquisition_date (datetime): The date of acquisition of the scenes.
    rel_orbit (int): The relative orbit of the scenes.

    Returns:
    tuple: A tuple containing the paths to the mosaicked B04 and cloud rasters.
    """

    # Create the extract path based on the base path, relative orbit, and acquisition date
    extract_path = os.path.join(base_path, f'R{rel_orbit:03}', f'{acquisition_date.strftime("%Y%m%d")}')

    # Define the output paths for the mosaicked B04 and cloud rasters
    b04_path_out = os.path.join(extract_path, f'S2-L2A-mosaic_{acquisition_date.strftime("%Y-%m-%dT%H%M%S")}_B04-10m')
    cld_path_out = os.path.join(extract_path, f'S2-L2A-mosaic_{acquisition_date.strftime("%Y-%m-%dT%H%M%S")}_cloud-20m')

    # Check if the B04 mosaic is missing
    if not os.path.exists(f'{b04_path_out}.vrt'):
        # Get the list of B04 files present in the extract path
        files_present_b04 = glob.glob(os.path.join(base_path, extract_path, '*_B04_10m.jp2'))

        # Create a list of files that are in UTM32
        files_present_b04_utm32 = []
        for file in files_present_b04:
            # Check if the file is in UTM32
            if utils.get_epsg(file) != 32632:
                # Warp the file to UTM32 using GDALWARP
                gsd_x, gsd_y = utils.get_pixel_spacing(file)
                os.system(
                    f'gdalwarp -tap -tr {gsd_x} {gsd_y} -q -ot UInt16 -s_srs EPSG:{utils.get_epsg(file)} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {file} {file.replace(".jp2", "_utm32.tif")}')
                files_present_b04_utm32.append(file.replace('.jp2', '_utm32.tif'))
            else:
                files_present_b04_utm32.append(file)

        # Create a VRT file for the B04 files
        print(f'\t- Creating vrt ({os.path.basename(b04_path_out)}.vrt)')
        os.system(
            f'gdalbuildvrt -allow_projection_difference -q -srcnodata 0 {b04_path_out}.vrt {" ".join(files_present_b04_utm32)}')

        ## Mosaic the B04 files using GDALWARP
        # print(f'\t- Mosaiking ({os.path.basename(b04_path_out)}.tif)')
        # pixel_spacing_x, pixel_spacing_y = utils.get_pixel_spacing(files_present_b04[0])
        # os.system(
        #     f'gdalwarp -tr {pixel_spacing_x} {pixel_spacing_y} -r near -tap -q -ot UInt16 -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {b04_path_out}.vrt {b04_path_out}.tif')

        ## Remove the intermediate UTM32 files
        # if len(glob.glob(os.path.join(base_path, extract_path, "*_utm32.tif"))) > 0:
        #    os.system(f'rm {" ".join(glob.glob(os.path.join(base_path, extract_path, "*_utm32.tif")))}')

    # Check if the cloud mosaic is missing
    if not os.path.exists(f'{cld_path_out}.vrt'):
        # Get the list of cloud files present in the extract path
        files_present_cld = glob.glob(os.path.join(base_path, extract_path, '*_MSK_CLDPRB_20m.jp2'))

        # Create a list of files that are in UTM32
        files_present_cld_utm32 = []
        for file in files_present_cld:
            # Check if the file is in UTM32
            if utils.get_epsg(file) != 32632:
                # Warp the file to UTM32 using GDALWARP
                os.system(
                    f'gdalwarp -q -ot Byte -s_srs EPSG:{utils.get_epsg(file)} -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS -co BIGTIFF=YES {file} {file.replace(".jp2", "_utm32.tif")}')
                files_present_cld_utm32.append(file.replace('.jp2', '_utm32.tif'))
            else:
                files_present_cld_utm32.append(file)

        # Create a VRT file for the cloud files
        print(f'\t- Creating vrt ({os.path.basename(cld_path_out)}.vrt)')
        os.system(
            f'gdalbuildvrt -q -allow_projection_difference -q -srcnodata 0 {cld_path_out}.vrt {" ".join(files_present_cld_utm32)}')

        ## Mosaic the cloud files using GDALWARP
        # print(f'\t- Mosaiking ({os.path.basename(cld_path_out)}.tif)')
        # pixel_spacing_x, pixel_spacing_y = utils.get_pixel_spacing(files_present_b04[0])
        # os.system(
        #     f'gdalwarp -tr {pixel_spacing_x} {pixel_spacing_y} -r near -tap -q -ot Byte -t_srs EPSG:32632 -overwrite -srcnodata 0 -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 -co NUM_THREADS=ALL_CPUS {cld_path_out}.vrt {cld_path_out}.tif')


    # Return the paths to the mosaicked B04 and cloud rasters
    return f'{b04_path_out}.vrt', f'{cld_path_out}.vrt'


def get_extent_and_dimensions(filepath):
    """
    This function gets the extent and dimensions of a raster file.

    Parameters:
    filepath (str): The path to the raster file.

    Returns:
    tuple: A tuple containing the minimum x, maximum x, minimum y, maximum y, x size, and y size of the raster.
    """
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


def resample_raster(src_filename, match_filename, dst_filename):
    """
    This function resamples a raster to match the resolution and extent of another raster.

    Parameters:
    src_filename (str): The path to the source raster file.
    match_filename (str): The path to the raster file to match.
    dst_filename (str): The path to the output raster file.
    """
    # Get information from the match file
    match_ds = gdal.Open(match_filename, gdalconst.GA_ReadOnly)
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize
    match_ds = None

    # Create a temporary file to change the resolution
    temp_filename = "/tmp/temp.tif"
    res = match_geotrans[1]
    os.system(
        f'gdalwarp -overwrite -q -tr {res} {res} -tap -co COMPRESS=DEFLATE -co PREDICTOR=2 {src_filename} {temp_filename}')

    # Now cut the temporary file to match the extent of the match file
    os.system(
        f'gdalwarp -overwrite -q -te {match_geotrans[0]} {match_geotrans[3] - high * res} {match_geotrans[0] + wide * res} {match_geotrans[3]} -dstnodata 0 -co COMPRESS=DEFLATE -co PREDICTOR=2 {temp_filename} {dst_filename}')

    # Remove the temporary file
    os.system(f'rm {temp_filename}')
    return dst_filename


def equalize_extents(im_reference, im_target):
    """
    This function equalizes the extents of two raster files.

    Parameters:
    im_reference (str): The path to the reference raster file.
    im_target (str): The path to the target raster file.

    Returns:
    str: The path to the output raster file.
    """
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
    gsd_x, gsd_y = utils.get_pixel_spacing(im_target)
    kwargs = {"format": "GTiff", "outputBounds": [minx_target, miny_target, maxx_target, maxy_target], "xRes": gsd_x,
              "yRes": gsd_y, "targetAlignedPixels": True}
    gdal.Warp(outputFile, im_reference, **kwargs)
    return outputFile


def shifts_to_tif(CRL, outfolder=None, mosaicing=False):
    """
    This function stores X- and Y-Offsets as geotiffs.

    Parameters:
    CRL (object): An object containing information about the coregistration.
    outfolder (str): The path to the output folder. If not provided, it defaults to the directory of the output file.

    Returns:
    None
    """
    from osgeo import gdal
    options = ['COMPRESS=DEFLATE', 'PREDICTOR=2', 'BIGTIFF=YES']

    # If outfolder is not provided, set it to the directory of the output file
    if outfolder is None:
        outfolder = os.path.dirname(CRL.path_out)

    # Create the output folder if it doesn't exist
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    # Check if there are no valid GCPs
    if len(CRL.coreg_info['GCPList']) == 0:
        # Open the input file
        with gdal.Open(CRL.im2shift.filePath) as ds:
            # Get the number of columns and rows
            cols = ds.RasterXSize
            rows = ds.RasterYSize

            # Create a zero-filled array for the X shifts
            x_shifts = np.full((rows, cols), 0)

            # Get the geotransform and projection from the input file
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()

            # Create a new geotiff file for the X shifts
            driver = gdal.GetDriverByName('GTiff')
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dx.tif')),
                                  cols, rows, 1, gdal.GDT_Byte, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing X Shifts')
            band.WriteArray(x_shifts)

            # Create a new geotiff file for the Y shifts
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dy.tif')),
                                  cols, rows, 1, gdal.GDT_Byte, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing Y Shifts')  # Because both are zero-matrices
            band.WriteArray(x_shifts)
            del x_shifts
    else:
        # Calculate the X shifts
        print('Calculating X Shifts')
        if mosaicing:
            CRL.tiepoint_grid.CoRegPoints_table.X_SHIFT_M = CRL.tiepoint_grid.CoRegPoints_table.X_SHIFT_M.astype(
                'float32')
            # Creating synthetic CRL with lowered gsd to minimize RAM requirement
            CRL_tmp, CRL_backup = utils.change_resolution_CRL(CRL, gsd_new=20)
            x_shifts = CRL_tmp.tiepoint_grid.to_interpolated_raster(metric='X_SHIFT_M', lowres_spacing=5)
            # Resampling the intermediary raster to the original resolution
            ## Create an interpolation function
            interp_func = interp2d(np.arange(x_shifts.shape[1]), np.arange(x_shifts.shape[0]), x_shifts, kind='cubic')

            # Evaluate the interpolation function at the new coordinates
            CRL_tmp = utils.change_resolution_CRL(CRL, backup=CRL_backup)
            x_new = np.linspace(0, x_shifts.shape[1], CRL_tmp.tiepoint_grid.shift.arr.shape[1])
            y_new = np.linspace(0, x_shifts.shape[0], CRL_tmp.tiepoint_grid.shift.arr.shape[0])
            x_shifts = np.int16(1000 * interp_func(x_new, y_new))

        else:
            x_shifts = CRL.tiepoint_grid.to_interpolated_raster(metric='X_SHIFT_M')

        # Open the input file
        with gdal.Open(CRL.im2shift.filePath) as ds:
            # Get the number of columns and rows
            cols = ds.RasterXSize
            rows = ds.RasterYSize

            # Get the geotransform and projection from the input file
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()

            # Create a new geotiff file for the X shifts
            driver = gdal.GetDriverByName('GTiff')
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dx.tif')),
                                  cols, rows, 1, gdal.GDT_Int16, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing X Shifts')
            band.WriteArray(x_shifts)
            outDs = None
            del x_shifts

        # Calculate the Y shifts
        print('Calculating Y Shifts')
        if mosaicing:
            CRL.tiepoint_grid.CoRegPoints_table.Y_SHIFT_M = CRL.tiepoint_grid.CoRegPoints_table.Y_SHIFT_M.astype(
                'float32')
            # Creating synthetic CRL with lowered gsd to minimize RAM requirement
            CRL_tmp, CRL_backup = utils.change_resolution_CRL(CRL, gsd_new=20)
            y_shifts = CRL_tmp.tiepoint_grid.to_interpolated_raster(metric='Y_SHIFT_M', lowres_spacing=5)

            # Resampling the intermediary raster to the original resolution
            ## Create an interpolation function
            interp_func = interp2d(np.arange(y_shifts.shape[1]), np.arange(y_shifts.shape[0]), y_shifts, kind='cubic')

            # Evaluate the interpolation function at the new coordinates
            CRL_tmp = utils.change_resolution_CRL(CRL, backup=CRL_backup)
            x_new = np.linspace(0, y_shifts.shape[1], CRL_tmp.tiepoint_grid.shift.arr.shape[1])
            y_new = np.linspace(0, y_shifts.shape[0], CRL_tmp.tiepoint_grid.shift.arr.shape[0])
            y_shifts = np.int16(1000 * interp_func(x_new, y_new))

        else:
            y_shifts = CRL.tiepoint_grid.to_interpolated_raster(metric='Y_SHIFT_M')

        # Open the input file
        with gdal.Open(CRL.im2shift.filePath) as ds:
            # Get the number of columns and rows
            cols = ds.RasterXSize
            rows = ds.RasterYSize

            # Get the geotransform and projection from the input file
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()

            # Create a new geotiff file for the Y shifts
            driver = gdal.GetDriverByName('GTiff')
            outDs = driver.Create(os.path.join(outfolder, os.path.basename(CRL.path_out).replace('.tif', '_dy.tif')),
                                  cols, rows, 1, gdal.GDT_Int16, options)
            outDs.SetGeoTransform(geotransform)
            outDs.SetProjection(projection)
            band = outDs.GetRasterBand(1)
            print('Writing Y Shifts')
            band.WriteArray(y_shifts)
            outDs = None
            del y_shifts

def coregister_S2(file_path_in, cld_mask_path_in, folder_path_out, cloud_threshold=65, mosaicing=False):
    """
    Coregister Sentinel-2 images using AROSICS.

    Args:
        file_path_in (str): Path to the input Sentinel-2 image file.
        cld_mask_path_in (str): Path to the input cloud mask file.
        folder_path_out (str): Path to the output folder.
        cloud_threshold (int, optional): Cloud threshold percentage (default: 65).
        mosaicing (bool, optional): Whether to perform mosaicing (default: False).
    """

    # Set maximum points for coregistration
    max_points = 5000

    # Set target image path
    im_target = file_path_in

    # Set output cloud mask path based on mosaicing flag
    if mosaicing:
        cld_mask_path_out = cld_mask_path_in.replace(".vrt", "_bin.tif")
    else:
        cld_mask_path_out = cld_mask_path_in.replace(".tif", "_bin.tif")

    # Get target image dimensions
    x_size, y_size = get_extent_and_dimensions(im_target)[-2:]

    # Calculate grid resolution based on maximum points and image dimensions
    grid_res = round((((x_size * y_size) / max_points) ** .5) / 5) * 5

    # Set output file name and folder based on mosaicing flag
    if mosaicing:
        out_name = os.path.basename(file_path_in).replace('B04', 'registration_swiss').replace('vrt', 'tif')
    else:
        out_name = os.path.basename(file_path_in).replace('B04', 'registration_swiss').replace('jp2', 'tif')
    out_folder = os.path.join(folder_path_out, f'TiePoint_GridRes_{grid_res}x{grid_res}px')

    # Check if output files already exist
    if not os.path.exists(os.path.join(out_folder, out_folder)) or not os.path.exists(
            os.path.join(out_folder, out_name.replace(".tif", ".shp"))) or not os.path.exists(
        os.path.join(out_folder, out_name.replace(".tif", "_dx.tif"))) or not os.path.exists(
        os.path.join(out_folder, out_name.replace(".tif", "_dy.tif"))):

        # Convert cloud probability to cloud mask
        if os.path.exists(cld_mask_path_out):
            gsd_b04 = utils.get_pixel_spacing(im_target)
            gsd_cld = utils.get_pixel_spacing(cld_mask_path_out)
            if not os.path.exists(cld_mask_path_out.replace(f'-{gsd_cld[0]:0.0f}m_', f'-{gsd_b04[0]:0.0f}m_')):
                cld_mask_path_out = resample_raster(cld_mask_path_out, im_target,
                                                    cld_mask_path_out.replace(f'-{gsd_cld[0]:0.0f}m_',
                                                                              f'-{gsd_b04[0]:0.0f}m_'))
            else:
                cld_mask_path_out = cld_mask_path_out.replace(f'-{gsd_cld[0]:0.0f}m_', f'-{gsd_b04[0]:0.0f}m_')
        else:
            print(f'\t- Converting Cloud Probability to Cloud Mask (≥{cloud_threshold}%)')
            os.system(
                f'gdal_calc.py -A {cld_mask_path_in} --overwrite --outfile={cld_mask_path_out} --calc="logical_and(A>={cloud_threshold}, A<=100)" --type=Byte --NoDataValue=0 --co COMPRESS=DEFLATE --co PREDICTOR=2 --co NUM_THREADS=ALL_CPUS --quiet')

            # Remove the intermediate UTM32 files
            if len(glob.glob(os.path.join(os.path.dirname(cld_mask_path_in), "*CLDPRB*_utm32.tif"))) > 0:
                os.system(
                    f'rm {" ".join(glob.glob(os.path.join(os.path.dirname(cld_mask_path_in), "*CLDPRB*_utm32.tif")))}')

            gsd_b04 = utils.get_pixel_spacing(im_target)
            gsd_cld = utils.get_pixel_spacing(cld_mask_path_out)
            cld_mask_path_out = resample_raster(cld_mask_path_out, im_target,
                                                cld_mask_path_out.replace(f'-{gsd_cld[0]:0.0f}m_',
                                                                          f'-{gsd_b04[0]:0.0f}m_'))

        # Set reference image path
        im_reference = '/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/S2_GRI.tif'

        # Print title
        title_str = f'{os.path.basename(im_target)} | {grid_res}x{grid_res}px'
        print('=' * len(title_str))
        print(title_str)
        print('-' * len(title_str))

        # Equalize extents of reference and target images
        im_reference_masked = equalize_extents(im_reference, im_target)

        # Limiting number of CPU threads to mitigate risk of out-of-memory error if using a mosaic
        if mosaicing:
            file_count = len(
                glob.glob(os.path.join(folder_path_out, 'S2*_MSIL2A_*_B04*10m.jp2')))  # Finding number of tiles
            num_cpus = utils.determine_cores(file_count)
            # num_cpus = 24 +
        else:
            num_cpus = max(os.cpu_count() - 1, 1)

        # Define coregistration arguments
        kwargs = {
            'grid_res': grid_res,
            'path_out': f'{out_name}',
            #'path_out': None, # Coregistered image not written
            'projectDir': out_folder,
            'q': False,
            'nodata': (0, 0),
            'mask_baddata_tgt': cld_mask_path_out,
            'out_crea_options': ['COMPRESS=DEFLATE', 'PREDICTOR=2'],
            'CPUs': num_cpus,
            'progress': False,
            'nodata': [0, 0],
        }

        # Perform coregistration
        CRL = COREG_LOCAL(im_reference_masked, im_target, **kwargs)
        if len(CRL.coreg_info['GCPList']) != 0:  # If there are GCPs, so not totally cloud covered
            CRL.correct_shifts()
            CRL.tiepoint_grid.to_PointShapefile(
                path_out=os.path.join(folder_path_out, f'TiePoint_GridRes_{grid_res}x{grid_res}px',
                                      out_name.replace(".tif", ".shp")))

        shifts_to_tif(CRL, mosaicing=mosaicing)

        print('=' * len(title_str))

    else:
        print('\t- Image already coregistered --> Skipping')


def mosaic_tiles(input_folder, operation, bands):
    """
    Mosaic Sentinel-2 tiles into a single image.

    Args:
        input_folder (str): Path to the input folder containing tile files.
        operation (str): Mosaicking operation to perform (e.g. 'average', 'median', etc.).
        bands (str): Bands to mosaic (e.g. 'both', 'dx', 'dy').
    """

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
        datestr_dt = datetime.datetime.strptime(dx_basename.split('_')[2], '%Y%m%dT%H%M%S')
        datestr = datestr_dt.strftime('%Y-%m-%dt%H%M%S')
    else:
        dy_basename = os.path.basename(dy_files[1])
        datestr_dt = datetime.datetime.strptime(dy_basename.split('_')[2], '%Y%m%dT%H%M%S')
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


def combine_geotiffs(output_folder, output_filename, input_geotiff_paths, operation):
    """
    Combine GeoTIFFs using a specified operation.

    Args:
        output_folder (str): Path to the output folder.
        output_filename (str): Name of the output file.
        input_geotiff_paths (list): List of input GeoTIFF file paths.
        operation (str): Operation to perform (e.g. 'average', 'median', etc.).
    """

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
    """
    Returns the corresponding NumPy data type for a given GDAL data type.

    :param gdal_data_type: The GDAL data type to convert.
    :return: The corresponding NumPy data type.
    """
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
    """
    Calculates the IDW (Inverse Distance Weighting) interpolation for a given shapefile folder and reference GeoTIFF.

    :param shapefile_folder: The folder containing the shapefiles.
    :param reference_geotiff: The reference GeoTIFF file.
    :param output_raster: The output raster file.
    :param power: The power to use for the IDW calculation (default=2).
    :param max_distance: The maximum distance to consider for the IDW calculation (default=None).
    :param max_points: The maximum number of points to consider for the IDW calculation (default=None).
    :param min_reliability: The minimum reliability required for a point to be considered (default=60).
    """
    import os
    import numpy as np
    from osgeo import ogr
    import rasterio, time

    # Open the reference GeoTIFF file
    with rasterio.open(reference_geotiff) as src:
        transform = src.transform
        rows, cols = src.height, src.width
        x_coords_px, y_coords_px = rasterio.transform.xy(transform, range(cols),
                                                         range(rows))  # Referenced to pixel center

    # Initialize arrays to store the x and y shift values
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

    # Filter features based on outlier and reliability criteria
    filtered_features = [feature for feature in merged_features if
                         feature.GetField("OUTLIER") != 1 and feature.GetField("RELIABILIT") >= min_reliability]
    del merged_features

    # Extract coordinates and shift values from filtered features
    x_coords_pts = np.array([feature.GetGeometryRef().GetX() for feature in filtered_features])
    y_coords_pts = np.array([feature.GetGeometryRef().GetY() for feature in filtered_features])
    x_shift_pts = np.array([feature.GetField('X_SHIFT_M') for feature in filtered_features])
    y_shift_pts = np.array([feature.GetField('Y_SHIFT_M') for feature in filtered_features])

    # Start timer
    start = time.time()

    # Iterate over x and y coordinates
    for x_count, x in enumerate(x_coords_px):
        if x_count > 0:
            print(
                f'{x_count / len(x_coords_px)} in {time.time() - start} sec (--> {((time.time() - start) / x_count) * len(x_coords_px)} remaining)')

        for y_count, y in enumerate(y_coords_px):
            # Calculate masks for x and y coordinates
            x_mask = np.abs(x - x_coords_pts) <= max_distance
            y_mask = np.abs(y - y_coords_pts) <= max_distance
            indices = np.where(x_mask & y_mask)[0]
            if len(indices) > 0:
                # Extract temporary arrays for IDW calculation
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
