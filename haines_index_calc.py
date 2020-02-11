#!/usr/bin/env python

# ------------------------------------------------------
# Calcolo dell'indice di Haines
# A lower atmosphere severity index for wildlife fires
#
# LOW ELEVATION: 950 - 850
# (A - lapse rate) = 950 mbT° - 850 mbT°
# (B - moisture) = 850 mbT° - 850 mbT° d
#
# MID ELEVATION: 850 - 700
# (A - lapse rate) = 850 mbT° - 700 mbT°
# (B - moisture) = 850 mbT° - 850 mbT° d
#
# HIGH ELEVATION: 700 - 500
# (A - lapse rate) = 700 mbT° - 500 mbT°
# (B - moisture) = 700 mbT° - 700 mbT° d
#
# ------------------------------------------------------

import sys
import getopt
import numpy as np
from osgeo import gdal
from datetime import datetime
import time
start_time = time.time()

driver = gdal.GetDriverByName("GTiff")
wkt_projection = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

def write_geotiff(filename, tempo, src_ds, array, type):
    """Scrivo i risultati in geotiff"""
    [rows, cols] = array.shape
    dataset = driver.Create(
        filename + '_' + type + '_' + tempo,
        cols,  # cols
        rows,  # rows
        1,
        gdal.GDT_Float64)

    dataset.SetGeoTransform(src_ds.GetGeoTransform())
    dataset.SetProjection(wkt_projection)

    dataset.GetRasterBand(1).WriteArray(array)

    dataset.FlushCache()
    dataset = None

def checktime(grib_time):
    """Estraggo la data per nominare i forecast"""
    data_hours = int(grib_time[:-8])
    return datetime.utcfromtimestamp(data_hours).strftime('%Y%m%d') + 'T' + datetime.utcfromtimestamp(
        data_hours).strftime('%H%M%S') + '000Z' + '.tiff'

def dewpoint_temp_calc(liv, temperature, umidity):
    """Calcolo la dew-point temperature per il livello"""
    pc = liv * 10000

    """definisco la rh al livello impostato"""
    exp = (0.7859 + 0.03477 * temperature) / (1.0 + 0.00412 * temperature) + 2
    rhc = pc * umidity / pow(10, exp) / (0.378 * umidity + 0.622)

    """definisco la td al livello impostato"""
    return temperature - ((14.55 + 0.114 * temperature) * (1 - 0.01 * rhc) + pow(
        (2.5 + 0.007 * temperature) * (1 - 0.01 * rhc), 3) + (
                                                  15.9 + 0.117 * temperature) * pow((1 - 0.01 * rhc),
                                                                                                     14))
def check_haines_index_type(type):
    if(type == 'low'):
        return [950, 850]
    if(type == 'mid'):
        return [850, 700]
    if(type == 'high'):
        return [700, 500]

def haines_index_calc(type):
    [sup, inf] = check_haines_index_type(type)
    src_ds = gdal.Open('../ecm.0p10.run00.grb')

    bands_count = src_ds.RasterCount

    temperature_sup_dataset_array = None
    temperature_inf_dataset_array = None
    specific_humidity_dataset_array = None

    for band in range(bands_count):
        band += 1
        src_subds = src_ds.GetRasterBand(band)
        metadata = src_subds.GetMetadata()

        if (metadata['GRIB_COMMENT'].find('Temperature [C]') != -1 and metadata['GRIB_SHORT_NAME'].find(str(sup)) != -1):
            temperature_sup_dataset = src_subds
            temperature_sup_dataset_array = temperature_sup_dataset.ReadAsArray()

        if (metadata['GRIB_COMMENT'].find('Temperature [C]') != -1 and metadata['GRIB_SHORT_NAME'].find(str(inf)) != -1):
            temperature_inf_dataset = src_subds
            temperature_inf_dataset_array = temperature_inf_dataset.ReadAsArray()

        """Estraggo la corretta specific humidity a seconda del tipo di elevation index"""
        if (type == 'low'):
            if (metadata['GRIB_COMMENT'].find('Specific humidity [kg/kg]') != -1 and metadata['GRIB_SHORT_NAME'].find(str(inf)) != -1):
                specific_humidity_dataset = src_subds
                specific_humidity_dataset_array = specific_humidity_dataset.ReadAsArray()
        else:
            if (metadata['GRIB_COMMENT'].find('Specific humidity [kg/kg]') != -1 and metadata['GRIB_SHORT_NAME'].find(str(sup)) != -1):
                specific_humidity_dataset = src_subds
                specific_humidity_dataset_array = specific_humidity_dataset.ReadAsArray()

        if (specific_humidity_dataset_array is not None and
                temperature_inf_dataset_array is not None and
                temperature_sup_dataset_array is not None):

            """Estraggo il tempo"""
            tempo = checktime(metadata['GRIB_VALID_TIME'])

            """Calcolo la dew-point temperatura per il livello a seconda del tipo di elevation index"""
            if (type == 'low'):
                tdc = dewpoint_temp_calc(inf, temperature_inf_dataset_array, specific_humidity_dataset_array)
            else:
                tdc = dewpoint_temp_calc(sup, temperature_sup_dataset_array, specific_humidity_dataset_array)


            """
            Calcolo 0000 GMT lapse rate
            (A) 950 mb T° - 850 mb T°
            """
            lapse_rate = temperature_sup_dataset_array - temperature_inf_dataset_array

            """Write geotiff lapse_rate_values"""
            write_geotiff('haines_images/lapse_rate_values', tempo, src_ds, lapse_rate, type)

            """Calcolo Factor Values (A)"""
            lapse_rate1_temp = np.less(lapse_rate, 4)
            np.putmask(lapse_rate, lapse_rate1_temp, 1)

            lapse_rate3_temp = np.greater_equal(lapse_rate, 8)
            np.putmask(lapse_rate, lapse_rate3_temp, 3)

            lapse_rate2_temp = np.logical_and([np.greater_equal(lapse_rate, 4)],
                                   [np.less(lapse_rate, 8)])
            np.putmask(lapse_rate, lapse_rate2_temp, 2)

            """Write geotiff lapse_rate_reclass"""
            write_geotiff('haines_images/lapse_rate_reclass', tempo, src_ds, lapse_rate, type)


            """
            Calcolo moisture a seconda del tipo di elevation index
            (B) 850 mb T° - 850 mb T° d
            """
            if (type == 'low'):
                moisture = temperature_inf_dataset_array - tdc
            else:
                moisture = temperature_sup_dataset_array - tdc

            """Write geotiff moisture_values"""
            write_geotiff('haines_images/moisture_values', tempo, src_ds, moisture, type)

            """Calcolo Factor Values (B)"""
            moisture1_temp = np.less(moisture, 6)
            np.putmask(moisture, moisture1_temp, 1)

            moisture3_temp = np.greater_equal(moisture, 10)
            np.putmask(moisture, moisture3_temp, 3)

            moisture2_temp = np.logical_and([np.greater_equal(moisture, 6)],
                                   [np.less(moisture, 10)])
            np.putmask(moisture, moisture2_temp, 2)

            """Write geotiff moisture_reclass"""
            write_geotiff('haines_images/moisture_reclass', tempo, src_ds, moisture, type)


            """
            Calcolo haines index
            factor values (A + B)
            """
            haines_index = lapse_rate + moisture

            """Write geotiff haines index values"""
            write_geotiff('haines_images/haines_index_values', tempo, src_ds, haines_index, type)

            """Class of day (potential for large fire)"""
            haines_index_verylow_temp = np.logical_or([np.equal(haines_index, 2)],
                                   [np.equal(haines_index, 3)])
            np.putmask(haines_index, haines_index_verylow_temp, 1)

            haines_index_low_temp = np.equal(haines_index, 4)
            np.putmask(haines_index, haines_index_low_temp, 2)

            haines_index_moderate_temp = np.equal(haines_index, 5)
            np.putmask(haines_index, haines_index_moderate_temp, 3)

            haines_index_high_temp = np.equal(haines_index, 6)
            np.putmask(haines_index, haines_index_high_temp, 4)

            """Write geotiff haines index"""
            write_geotiff('haines_images/haines_index_reclass', tempo, src_ds, haines_index, type)


            temperature_sup_dataset_array = None
            temperature_inf_dataset_array = None
            specific_humidity_dataset_array = None

def print_usage():
    print("haines_index_calc.py -e <elevation>")

def main(argv):
    #print("ARGV      :", sys.argv[1:])
    try:
        # opts is a list of returning key-value pairs, args is the options left after striped
        # the short options 'hi:o:', if an option requires an input, it should be followed by a ":"
        # the long options 'ifile=' is an option that requires an input, followed by a "="
        opts, args = getopt.getopt(argv, "he:", ["help", "elevation="])
    except getopt.GetoptError as err:
        print(err)  # will print something like "option -a not recognized"
        print_usage()
        sys.exit(2)
    if not opts:
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("haines_index_calc.py -e <elevation>")
            sys.exit(2)
        elif opt in ("-e", "--elevation"):
            elev = arg
        else:
            assert False, "unhandled option"

    try:
        haines_index_calc(elev)
    except Exception as e:
        print(e)
        sys.exit(2)

if __name__ == '__main__':
    main(sys.argv[1:])
    print("--- %s seconds ---" % (time.time() - start_time))