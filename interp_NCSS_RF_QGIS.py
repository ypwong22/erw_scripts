from qgis.core import QgsRasterLayer, QgsProject
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
import processing
import os

for layer in range(10):

    input_raster_path = os.path.join(os.getcwd(), f'RF_NCSS_layer_{layer}.tif')
    output_path = os.path.join(os.getcwd())

    band_count = 6

    ###########################################################################
    # Process each band
    ###########################################################################
    filled_bands = []
    filled_ucs = []

    input_raster = QgsRasterLayer(input_raster_path, f"RF_NCSS_layer_{layer}")
    QgsProject.instance().addMapLayer(input_raster)

    for band in range(1, band_count + 1):
        print(f"Processing band {band}")

        # Select a band from the input_raster_path
        output_single_band = os.path.join(os.getcwd(), f'temp_input_{band}.tif')

        raster_entries = []
        raster_entry = QgsRasterCalculatorEntry()
        raster_entry.raster = input_raster
        raster_entry.bandNumber = band
        raster_entry.ref = 'input@' + str(band)
        raster_entries.append(raster_entry)

        calc = QgsRasterCalculator(f'{raster_entry.ref}', 
                                    output_single_band, 
                                    'GTiff', 
                                    input_raster.extent(), 
                                    input_raster.width(), 
                                    input_raster.height(), 
                                    raster_entries)
        calc.processCalculation()

        single_band_raster = QgsRasterLayer(output_single_band, 'single_band')

        # Fill no data
        filled_output = os.path.join(output_path, f"filled_band_{band}.tif")
        filled_uc = os.path.join(output_path, f"filled_uc_{band}.tif")
        processing.run("grass7:r.fill.stats", {
            'input': single_band_raster,
            '-k': True,
            'mode': 0,
            '-m': False,
            'distance': 3,
            'minimum': None,
            'maximum': None,
            'power': 2,
            'cells': 8,
            'output': filled_output,
            'uncertainty': filled_uc, 
            'GRASS_REGION_PARAMETER': None,
            'GRASS_REGION_CELLSIZE_PARAMETER': 0,
            'GRASS_RASTER_FORMAT_OPT': '',
            'GRASS_RASTER_FORMAT_META': ''
        })

        filled_bands.append(filled_output)
        filled_ucs.append(filled_uc)

        single_band_raster = None
        filled_output = None
        filled_uc = None

    ###########################################################################
    # Merge bands into a single multi-band raster
    ###########################################################################
    processing.run("gdal:merge", {
        'INPUT': filled_bands,
        'PCT': False,
        'SEPARATE': True,
        'NODATA_INPUT': -9999,
        'NODATA_OUTPUT': -9999,
        'OPTIONS': '',
        'EXTRA': '',
        'DATA_TYPE': 5,
        'OUTPUT': os.path.join(output_path, f'RF_NCSS_layer_{layer}_filled.tif')
    })

    processing.run("gdal:merge", {
        'INPUT': filled_ucs,
        'PCT': False,
        'SEPARATE': True,
        'NODATA_INPUT': -9999,
        'NODATA_OUTPUT': -9999,
        'OPTIONS': '',
        'EXTRA': '',
        'DATA_TYPE': 5,
        'OUTPUT': os.path.join(output_path, f'RF_NCSS_layer_{layer}_uc.tif')
    })

    output_raster = QgsRasterLayer(os.path.join(output_path, f'RF_NCSS_layer_{layer}_filled.tif'), 
                                f"RF_NCSS_layer_{layer}_filled")
    QgsProject.instance().addMapLayer(output_raster)

    print(f"Processing complete. Output saved to: {output_path}")

for band in range(1, band_count+1):
    os.remove(os.path.join(output_path, f'temp_input_{band}.tif'))
    os.remove(os.path.join(output_path, f'temp_input_{band}.tif.aux.xml'))
    os.remove(os.path.join(output_path, f'filled_band_{band}.tif'))
    os.remove(os.path.join(output_path, f'filled_band_{band}.tfw'))
    os.remove(os.path.join(output_path, f'filled_band_{band}.tif.aux.xml'))
    os.remove(os.path.join(output_path, f'filled_uc_{band}.tif'))
    os.remove(os.path.join(output_path, f'filled_uc_{band}.tfw'))
    os.remove(os.path.join(output_path, f'filled_uc_{band}.tif.aux.xml'))
