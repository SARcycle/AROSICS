import glob
import os
import subprocess

import numpy as np
from osgeo import gdal


def create_timelapse_gif(image_path, extent, output_folder, output_filename):
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Filter images that have at least one non-null pixel within the extent
    images_to_stack = []
    for file in glob.glob(image_path):
        ds = gdal.Open(file)
        # Get the geotransform matrix
        gt = ds.GetGeoTransform()

        # Get the extent of the image
        x_min_tif = gt[0]
        y_max_tif = gt[3]
        x_max_tif = x_min_tif + gt[1] * ds.RasterXSize
        y_min_tif = y_max_tif + gt[5] * ds.RasterYSize

        x_min_extent, y_min_extent, x_max_extent, y_max_extent = np.round(np.array(extent) / gt[1]) * gt[1]

        x_axis_tif = np.arange(start=x_min_tif, stop=x_max_tif, step=gt[1])
        y_axis_tif = np.arange(start=y_max_tif, stop=y_min_tif, step=gt[5])
        x_min_idx = np.argmin(np.abs(x_axis_tif - x_min_extent))
        x_max_idx = np.argmin(np.abs(x_axis_tif - x_max_extent))
        y_min_idx = np.argmin(np.abs(y_axis_tif - y_min_extent))
        y_max_idx = np.argmin(np.abs(y_axis_tif - y_max_extent))

        cropped_array = ds.GetRasterBand(1).ReadAsArray(x_min_idx, y_max_idx, x_max_idx - x_min_idx,
                                                        y_min_idx - y_max_idx)
        no_data_value = ds.GetRasterBand(1).GetNoDataValue()
        if np.any((cropped_array != 0) & (cropped_array != no_data_value)):
            images_to_stack.append(file)

    # Sort the images by date (assuming date is in the penultimate subfolder in the format YYYYMMDD)
    images_to_stack.sort(key=lambda x: os.path.basename(os.path.dirname(os.path.dirname(x))))

    # Resize the images to 1080p and create a timelapse gif using FFmpeg
    resized_images = []
    for i, file in enumerate(images_to_stack):
        ds = gdal.Open(file)
        # Get the geotransform matrix
        gt = ds.GetGeoTransform()

        # Get the extent of the image
        x_min_tif = gt[0]
        y_max_tif = gt[3]
        x_max_tif = x_min_tif + gt[1] * ds.RasterXSize
        y_min_tif = y_max_tif + gt[5] * ds.RasterYSize

        x_min_extent, y_min_extent, x_max_extent, y_max_extent = np.round(np.array(extent) / gt[1]) * gt[1]

        x_axis_tif = np.arange(start=x_min_tif, stop=x_max_tif, step=gt[1])
        y_axis_tif = np.arange(start=y_max_tif, stop=y_min_tif, step=gt[5])
        x_min_idx = np.argmin(np.abs(x_axis_tif - x_min_extent))
        x_max_idx = np.argmin(np.abs(x_axis_tif - x_max_extent))
        y_min_idx = np.argmin(np.abs(y_axis_tif - y_min_extent))
        y_max_idx = np.argmin(np.abs(y_axis_tif - y_max_extent))

        output_file = f"{i:03d}_resized.tif"
        # driver = gdal.GetDriverByName('GTiff')
        # dst_ds = driver.CreateCopy(output_file, ds, 0, [x_min_idx, y_max_idx, x_max_idx - x_min_idx, y_min_idx - y_max_idx])
        # dst_ds = None
        convert_geotiff(file, f"{i:03d}_resized_1080p.tif",
                        (x_min_idx, y_max_idx, x_max_idx - x_min_idx, y_min_idx - y_max_idx))
        resized_images.append(f"{i:03d}_resized_1080p.tif")

    # Create the timelapse gif
    print(f'Creating timelapse {output_folder}/{output_filename}')
    subprocess.run(
        ['ffmpeg', '-framerate', '2', '-i', '%03d_resized_1080p.tif', '-y', f'{output_folder}/{output_filename}'])

    # Clean up the temporary files
    for file in resized_images:
        os.remove(file)


def convert_geotiff(input_file, output_file, crop_extent):
    import rasterio
    import numpy as np
    from PIL import Image
    print(f'Resizing {input_file}')
    with rasterio.open(input_file) as src:
        # Read the image data within the crop extent
        img = src.read(1, window=rasterio.windows.Window(*crop_extent))

        # Calculate the 2% and 98% percentiles of the histogram
        p2, p98 = np.percentile(img, (2, 98))

        # Clip the image values to the calculated percentiles
        clipped_img = np.clip(img, p2, p98)

        # Scale the clipped image to the range [0, 255] for JPEG conversion
        scaled_img = (clipped_img - p2) / (p98 - p2) * 255
        scaled_img = scaled_img.astype(np.uint8)

        # Rescale the image to a maximum of 1920 pixels in width or 1080 pixels in height
        max_width, max_height = 1920, 1080
        width, height = scaled_img.shape[1], scaled_img.shape[0]
        scale_factor = min(max_width / width, max_height / height)
        new_width, new_height = int(width * scale_factor), int(height * scale_factor)

        # Use Pillow to rescale the image with better interpolation
        img_pil = Image.fromarray(scaled_img)
        rescaled_img = img_pil.resize((new_width, new_height), Image.NEAREST)

        # Save the clipped and scaled image as a JPEG
        rescaled_img.save(output_file)

        return output_file


# Example usage:
extent_schoenbuehl = ((384045, 5207713, 386346, 5208780), 'Schoenbuehl')
extent_rapperswil = ((484019, 5228534, 487054, 5230748), 'Rapperswil')
extent_chur = ((537250, 5188630, 540285, 5190872), 'Chur')
extent_poschiavo = ((583372, 5123844, 586423, 5126108), 'Poschiavo')
extent_martigny = ((349936, 5107514, 353078, 5109648), 'Martigny')
extent_geneve = ((274551, 5121786, 277722, 5123897), 'Geneve')

orbit = 'R*'

for orbit in ['R008', 'R022', 'R065', 'R108', 'R*']:

    date = '*'

    if orbit == 'R*':
        orbit_name = 'allOrbits'
    else:
        orbit_name = orbit
    if date == '*':
        date_name = 'allDates'
    elif '*' in date:
        date_name = date.replace('*', '')
    else:
        date_name = date

    extent = extent_schoenbuehl

    image_path = f'/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2/{orbit}/{date}/TiePoint_GridRes_*x*px/S2-L2A-mosaic_*_registration_swiss-10m.tif'
    output_folder = '/home/localadmin/Downloads'
    output_filename = f'timelapse_{extent[1]}_{orbit_name}_{date_name}.mp4'
    create_timelapse_gif(image_path, extent[0], output_folder, output_filename)

    extent = extent_rapperswil

    image_path = f'/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2/{orbit}/{date}/TiePoint_GridRes_*x*px/S2-L2A-mosaic_*_registration_swiss-10m.tif'
    output_folder = '/home/localadmin/Downloads'
    output_filename = f'timelapse_{extent[1]}_{orbit_name}_{date_name}.mp4'
    create_timelapse_gif(image_path, extent[0], output_folder, output_filename)

    extent = extent_chur

    image_path = f'/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2/{orbit}/{date}/TiePoint_GridRes_*x*px/S2-L2A-mosaic_*_registration_swiss-10m.tif'
    output_folder = '/home/localadmin/Downloads'
    output_filename = f'timelapse_{extent[1]}_{orbit_name}_{date_name}.mp4'
    create_timelapse_gif(image_path, extent[0], output_folder, output_filename)

    extent = extent_poschiavo

    image_path = f'/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2/{orbit}/{date}/TiePoint_GridRes_*x*px/S2-L2A-mosaic_*_registration_swiss-10m.tif'
    output_folder = '/home/localadmin/Downloads'
    output_filename = f'timelapse_{extent[1]}_{orbit_name}_{date_name}.mp4'
    create_timelapse_gif(image_path, extent[0], output_folder, output_filename)

    extent = extent_martigny

    image_path = f'/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2/{orbit}/{date}/TiePoint_GridRes_*x*px/S2-L2A-mosaic_*_registration_swiss-10m.tif'
    output_folder = '/home/localadmin/Downloads'
    output_filename = f'timelapse_{extent[1]}_{orbit_name}_{date_name}.mp4'
    create_timelapse_gif(image_path, extent[0], output_folder, output_filename)

    extent = extent_geneve

    image_path = f'/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2/{orbit}/{date}/TiePoint_GridRes_*x*px/S2-L2A-mosaic_*_registration_swiss-10m.tif'
    output_folder = '/home/localadmin/Downloads'
    output_filename = f'timelapse_{extent[1]}_{orbit_name}_{date_name}.mp4'
    create_timelapse_gif(image_path, extent[0], output_folder, output_filename)
