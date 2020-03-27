#!/usr/bin/env python

# ------------------------------------------------------
# Calcolo dell'indice di Haines
# A lower atmosphere severity index for wildlife fires
#
# LOW ELEVATION: 950 - 850 PER QUOTE < 300 m
# (A - lapse rate) = 950 mbT° - 850 mbT°
# (B - moisture) = 850 mbT° - 850 mbT° d
#
# MID ELEVATION: 850 - 700 PER QUOTE >= 300 m e QUOTE < 900
# (A - lapse rate) = 850 mbT° - 700 mbT°
# (B - moisture) = 850 mbT° - 850 mbT° d
#
# HIGH ELEVATION: 700 - 500 PER QUOTE >= 900
# (A - lapse rate) = 700 mbT° - 500 mbT°
# (B - moisture) = 700 mbT° - 700 mbT° d
#
# ------------------------------------------------------

"""
https://www.earthdatascience.org/courses/earth-analytics-bootcamp/functions/apply-functions-numpy-arrays/
"""

import sys
import getopt
import numpy as np
from osgeo import gdal
from datetime import datetime
import time
start_time = time.time()

driver = gdal.GetDriverByName("GTiff")
wkt_projection = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

def write_geotiff(filename, tempo, src_ds, array, type_str, type = gdal.GDT_Float64):
    """Scrivo i risultati in geotiff"""
    [rows, cols] = array.shape
    dataset = driver.Create(
        filename + '_' + type_str + '_' + tempo + '.tiff',
        cols,  # cols
        rows,  # rows
        1,
        type)

    dataset.SetGeoTransform(src_ds.GetGeoTransform())
    dataset.SetProjection(wkt_projection)

    dataset.GetRasterBand(1).WriteArray(array)

    dataset.FlushCache()
    dataset = None

def checktime(grib_time):
    """Estraggo la data per nominare i forecast"""
    data_hours = int(grib_time[:-8])
    return datetime.utcfromtimestamp(data_hours).strftime('%Y%m%d') + 'T' + datetime.utcfromtimestamp(
        data_hours).strftime('%H%M%S') + '000Z'

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

def write_haines_array(key, value, src_ds, variable_dict, i, writegeotiff, variable_haines_list):

        """Estraggo il tempo"""
        tempo = variable_dict['tempi'][i]

        """
        (A) Calcolo 0000 GMT lapse rate
        """
        lapse_rate = variable_dict['temperature_sup_dataset_array_dict'][key][tempo] - variable_dict['temperature_inf_dataset_array_dict'][key][tempo]

        """Write geotiff lapse_rate_values"""
        if (writegeotiff):
            write_geotiff('haines_images/lapse_rate_values', tempo, src_ds, lapse_rate, key)

        """Calcolo Factor Values (A)"""
        lapse_rate1_temp = np.less(lapse_rate, 4)
        np.putmask(lapse_rate, lapse_rate1_temp, 1)

        lapse_rate3_temp = np.greater_equal(lapse_rate, 8)
        np.putmask(lapse_rate, lapse_rate3_temp, 3)

        lapse_rate2_temp = np.logical_and([np.greater_equal(lapse_rate, 4)],
                                          [np.less(lapse_rate, 8)])
        np.putmask(lapse_rate, lapse_rate2_temp, 2)

        """Write geotiff lapse_rate_reclass"""
        if (writegeotiff):
            write_geotiff('haines_images/lapse_rate_reclass', tempo, src_ds, lapse_rate, key)

        """Calcolo la dew-point temperatura per il livello a seconda del tipo di elevation index"""

        """If LOW calcolo la Temperatura di rugiada del livello più BASSO"""
        tdc_low = dewpoint_temp_calc(value['inf'], variable_dict['temperature_inf_dataset_array_dict'][key][tempo], variable_dict['specific_humidity_inf_dataset_array_dict'][key][tempo])

        """If MID or HIGH calcolo la Temperatura di rugiada del livello più ALTO"""
        tdc_mid_high = dewpoint_temp_calc(value['sup'], variable_dict['temperature_sup_dataset_array_dict'][key][tempo], variable_dict['specific_humidity_sup_dataset_array_dict'][key][tempo])

        """
        (B) Calcolo moisture a seconda del tipo di elevation index
        """
        if (key == 'low'):
            moisture = variable_dict['temperature_inf_dataset_array_dict'][key][tempo] - tdc_low
        else:
            moisture = variable_dict['temperature_sup_dataset_array_dict'][key][tempo] - tdc_mid_high

        """Write geotiff moisture_values"""
        if (writegeotiff):
            write_geotiff('haines_images/moisture_values', tempo, src_ds, moisture, key)

        """Calcolo Factor Values (B)"""
        moisture1_temp = np.less(moisture, 6)
        np.putmask(moisture, moisture1_temp, 1)

        moisture3_temp = np.greater_equal(moisture, 10)
        np.putmask(moisture, moisture3_temp, 3)

        moisture2_temp = np.logical_and([np.greater_equal(moisture, 6)],
                                        [np.less(moisture, 10)])
        np.putmask(moisture, moisture2_temp, 2)

        """Write geotiff moisture_reclass"""
        if (writegeotiff):
            write_geotiff('haines_images/moisture_reclass', tempo, src_ds, moisture, key)

        """
        Calcolo haines index
        factor values (A + B)
        """
        haines_index = lapse_rate + moisture

        """Write geotiff haines index values"""
        if (writegeotiff):
            write_geotiff('haines_images/haines_index_values', tempo, src_ds, haines_index, key)

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
        write_geotiff('haines_images/haines_index_reclass', tempo, src_ds, haines_index, key, gdal.GDT_Int16)

        variable_haines_list[i]['haines'][key] = haines_index

        return variable_haines_list

def read_variables(key, value, src_subds, metadata, src_ds, variable_dict):
    tempo = checktime(metadata['GRIB_VALID_TIME'])
    if (metadata['GRIB_COMMENT'].find('Geopotential (at the surface = orography) [m^2/s^2]') != -1 and metadata[
        'GRIB_SHORT_NAME'].find('0-SFC') != -1):
        geopotential = src_subds
        geopotential_array = geopotential.ReadAsArray()
        variable_dict['geopotential_array_dict'] = geopotential_array / 9.80665
        write_geotiff('haines_images/orography', tempo, src_ds, geopotential_array / 9.80665, 'orography')

    if (metadata['GRIB_COMMENT'].find('Temperature [C]') != -1 and metadata['GRIB_SHORT_NAME'].find(str(value['sup'])) != -1):
        temperature_sup_dataset = src_subds
        temperature_sup_dataset_array = temperature_sup_dataset.ReadAsArray()
        variable_dict['temperature_sup_dataset_array_dict'][key][tempo] = temperature_sup_dataset_array

    if (metadata['GRIB_COMMENT'].find('Temperature [C]') != -1 and metadata['GRIB_SHORT_NAME'].find(str(value['inf'])) != -1):
        temperature_inf_dataset = src_subds
        temperature_inf_dataset_array = temperature_inf_dataset.ReadAsArray()
        variable_dict['temperature_inf_dataset_array_dict'][key][tempo] = temperature_inf_dataset_array

    """Estraggo la corretta specific humidity a seconda del tipo di elevation index"""
    if (metadata['GRIB_COMMENT'].find('Specific humidity [kg/kg]') != -1 and metadata['GRIB_SHORT_NAME'].find(
            str(value['inf'])) != -1):
        specific_humidity_inf_dataset = src_subds
        specific_humidity_inf_dataset_array = specific_humidity_inf_dataset.ReadAsArray()
        variable_dict['specific_humidity_inf_dataset_array_dict'][key][tempo] = specific_humidity_inf_dataset_array

    if (metadata['GRIB_COMMENT'].find('Specific humidity [kg/kg]') != -1 and metadata['GRIB_SHORT_NAME'].find(
            str(value['sup'])) != -1):
        specific_humidity_sup_dataset = src_subds
        specific_humidity_sup_dataset_array = specific_humidity_sup_dataset.ReadAsArray()
        variable_dict['specific_humidity_sup_dataset_array_dict'][key][tempo] = specific_humidity_sup_dataset_array
        variable_dict['tempi'].append(tempo)

    return variable_dict

def haines_index_calc(types):

    src_ds = gdal.Open('../ecm.0p10.run00.grb')

    bands_count = src_ds.RasterCount

    variable_dict = {
        'geopotential_array_dict': {},
        'temperature_sup_dataset_array_dict': {
            'low': {},
            'mid': {},
            'high': {}
        },
        'temperature_inf_dataset_array_dict': {
            'low': {},
            'mid': {},
            'high': {}
        },
        'specific_humidity_inf_dataset_array_dict': {
            'low': {},
            'mid': {},
            'high': {}
        },
        'specific_humidity_sup_dataset_array_dict': {
            'low': {},
            'mid': {},
            'high': {}
        },
        'tempi': []
    }

    for key, value in types.items():
        variable_dict['tempi'] = []
        for band in range(bands_count):
            band += 1
            src_subds = src_ds.GetRasterBand(band)
            metadata = src_subds.GetMetadata()

            variable_dict = read_variables(key, value, src_subds, metadata, src_ds, variable_dict)

    variable_dict['tempi'].sort()
    dict_len = len(variable_dict['tempi'])

    variable_haines_list = []

    for i in range(dict_len):

        variable_haines_list.append({
            'tempo': variable_dict['tempi'][i],
            'elev': variable_dict['geopotential_array_dict'],
            'haines': {
                'low': {},
                'mid': {},
                'high': {}
            }
        })

        for key, value in types.items():
            write_haines_array(key, value, src_ds, variable_dict, i, False, variable_haines_list)

    variable_haines_list_lenght = len(variable_haines_list)

    for c in range(variable_haines_list_lenght):
        geopotential_array_300 = variable_haines_list[c]['elev'].copy()
        geopotential_array_300_900 = variable_haines_list[c]['elev'].copy()
        geopotential_array_900 = variable_haines_list[c]['elev'].copy()


        geopotential_array_300_temp1 = np.less_equal(geopotential_array_300, 300)
        np.putmask(geopotential_array_300, geopotential_array_300_temp1, 1)

        geopotential_array_300_temp2 = np.greater(geopotential_array_300, 300)
        np.putmask(geopotential_array_300, geopotential_array_300_temp2, 0)
        #write_geotiff('haines_images_all/geopotential_array_300', variable_haines_list[c]['tempo'], src_ds, geopotential_array_300, 'LOW')


        geopotential_array_300_900_temp1 = np.less_equal(geopotential_array_300_900, 300)
        np.putmask(geopotential_array_300_900, geopotential_array_300_900_temp1, 0)

        geopotential_array_300_900_temp3 = np.logical_and([np.greater(geopotential_array_300_900, 300)],
                                        [np.less_equal(geopotential_array_300_900, 900)])
        np.putmask(geopotential_array_300_900, geopotential_array_300_900_temp3, 1)

        geopotential_array_300_900_temp2 = np.greater(geopotential_array_300_900, 900)
        np.putmask(geopotential_array_300_900, geopotential_array_300_900_temp2, 0)
        #write_geotiff('haines_images_all/geopotential_array_300_900', variable_haines_list[c]['tempo'], src_ds, geopotential_array_300_900, 'MID')


        geopotential_array_900_temp1 = np.less_equal(geopotential_array_900, 900)
        np.putmask(geopotential_array_900, geopotential_array_900_temp1, 0)

        geopotential_array_900_temp2 = np.greater(geopotential_array_900, 900)
        np.putmask(geopotential_array_900, geopotential_array_900_temp2, 1)
        #write_geotiff('haines_images_all/geopotential_array_900', variable_haines_list[c]['tempo'], src_ds, geopotential_array_900, 'HIGH')


        new_low = np.multiply(variable_haines_list[c]['haines']['low'], geopotential_array_300)
        #write_geotiff('haines_images_all/haines_index_reclass', variable_haines_list[c]['tempo'], src_ds, new_low, 'LOW')

        new_mid = np.multiply(variable_haines_list[c]['haines']['mid'], geopotential_array_300_900)
        #write_geotiff('haines_images_all/haines_index_reclass', variable_haines_list[c]['tempo'], src_ds, new_mid, 'MID')

        new_high = np.multiply(variable_haines_list[c]['haines']['high'], geopotential_array_900)
        #write_geotiff('haines_images_all/haines_index_reclass', variable_haines_list[c]['tempo'], src_ds, new_high,'HIGH')

        new_total = new_low+new_mid+new_high

        write_geotiff('haines_images_all/haines_index_reclass', variable_haines_list[c]['tempo'], src_ds, new_total, 'ALL')





"""
OROGRAFIA
The geopotential height can be calculated by dividing the geopotential by the Earth's gravitational acceleration, g (=9.80665 m s-2).
The geopotential height plays an important role in synoptic meteorology (analysis of weather patterns).
Charts of geopotential height plotted at constant pressure levels (e.g., 300, 500 or 850 hPa)
can be used to identify weather systems such as cyclones, anticyclones, troughs and ridges.
  Description = 0[-] SFC (Ground or water surface)
  Metadata:
    GRIB_COMMENT=Geopotential (at the surface = orography) [m^2/s^2]
    GRIB_ELEMENT=Z
    GRIB_FORECAST_SECONDS=75600 sec
    GRIB_REF_TIME=  1581292800 sec UTC
    GRIB_SHORT_NAME=0-SFC
    GRIB_UNIT=[m^2/s^2]
"""
def print_usage():
    print("haines_index_calc_all.py -e <elevation>")

def main(argv):
    #print("ARGV      :", sys.argv[1:])
    elev = None
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
            print("haines_index_calc_all.py -e <elevation>")
            sys.exit(2)
        elif opt in ("-e", "--elevation"):
            elev = arg
        else:
            assert False, "unhandled option"

    try:
        elev = {
            'low': {
                'sup': 950,
                'inf': 850
            },
            'mid': {
                'sup': 850,
                'inf': 700
            },
            'high': {
                'sup': 700,
                'inf': 500
            }
        }
        haines_index_calc(elev)
    except Exception as e:
        print(e)
        sys.exit(2)

if __name__ == '__main__':
    main(sys.argv[1:])
    print("--- %s seconds ---" % (time.time() - start_time))