#!/usr/bin/env python

"""
Calcolo della maschera finale VERO-FALSO per individuare
le finestre ambientali per l'applicazione del fuoco prescritto.

ctl file
dset ^/home/lamma/data/grib/all_ecm_3km/incendi/incendi_arw_ecm_3km_run00.gra
undef 9.999E+20
title Model Data for Fire
xdef 267 LINEAR 8 0.03
ydef 201 LINEAR 40 0.03
zdef 1 levels 1000
tdef 3 linear 00Z18FEB2020 1dy
vars 4
tmpsfc 0 11,1,0 ** Mean Daily T2m [C]
apcpsfc 0 61,1,0 ** Total Daily precipitation [kg/m^2]
rhsfc 0 52,1,0 ** Mean Daily 2m Relative Humidity [%]
wind10m 0 32,105,10 ** Mean Daily 10 m Wind Velocity [m/s]
endvars
"""

import os
import getopt
import sys
from sys import path

path.remove('/mnt/hd/sviluppo/library/libimage/gdal/swig/python/build/lib.linux-x86_64-2.7')

import numpy as np
from osgeo import gdal, gdalconst
from datetime import timedelta, datetime
import time
start_time = time.time()
import fnmatch
from functools import reduce
import subprocess
import logging

wkt_projection = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
wkt_projection_prec = 'PROJCS["WGS 84 / UTM zone 32N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32632"]]'
tmp_directory = '/mnt/hd/operativo/tmp/fuoco_prescr_tmp_data/'
driver_list = ["GTiff", "RST"]
driver = gdal.GetDriverByName(driver_list[1])
driver_ext = '.rst'

"""
    SOGLIE ORIGINALI

    "threshold_tmpsfc": [0, 18],
    "threshold_rhsfc": [40, 75],
    "threshold_wind10m": [1, 15],
    "threshold_ffmc": [80, 95],
    "threshold_dmc": [0, 20],
    "threshold_prec": [1, 7],
"""

threshold_dict = {
    "threshold_tmpsfc": [-5, 20],
    "threshold_rhsfc": [40, 80],
    "threshold_wind10m": [0, 15],
    "threshold_ffmc": [75, 95],
    "threshold_dmc": [0, 20],
    "threshold_prec": [3, 20],
}

def get_threshold(params):
    dict = threshold_dict
    threshold = list(dict.get(params))
    return threshold

"""
SetGeoTransform gdal parameters

1 - The Upper Left easting coordinate (i.e., horizontal)
2 - The E-W pixel spacing
3 - The rotation (0 degrees if image is "North Up")
4 - The Upper left northing coordinate (i.e., vertical)
5 - The rotation (0 degrees)
6 - The N-S pixel spacing, negative as we will be counting from the UL corner
"""
def get_geotransform(params):
    dict = {
        "model": {
            1: 8.0000000000000000,
            2: 0.02999999999999999542,
            3: 0,
            4: 40.0000000000000000,
            5: 0,
            6: 0.02999999999999999542
        },
        "rst": {
            1: 7.9850000000000003,
            2: 0.03005617977528090776,
            3: 0,
            4: 46.0299989999999966,
            5: 0,
            6: -0.03007461691542288526
        },
        "prec": {
            1: 548174.9000000000232831,
            2: 1000.0000000,
            3: 0,
            4: 4933721.0000000000000000,
            5: 0,
            6: -1000.0000000
        }
    }
    geotransform_list = list(dict.get(params).values())
    return geotransform_list

nlat = 201
nlon = 267
nz = 1

dt = np.dtype((np.float32, (nlat, nlon)))

def read(f, n=1, offset=0):
    """in numpy > 1.17.3"""
    return np.fromfile(f, dt, n, "", offset)

def calc_byte(n, var):
    img_bytes = nlon * nlat * 4
    if (n == 0):
        return img_bytes * var

    if (n == 1):
        return img_bytes * (var + 4)

    if (n == 2):
        return img_bytes * (var + 8)

    if (n == 3):
        return img_bytes * (var + 12)


def write_geotiff_file(data, filename, transformparams):
    if transformparams == 'model':
        """in numpy > 1.17.3"""
        data = np.squeeze(data, axis=0)
        if driver.ShortName == 'RST':
            """in numpy > 1.17.3"""
            data = np.flip(data, axis=0)

    [rows, cols] = data.shape

    dataset = driver.Create(
        filename,
        cols,
        rows,
        1,
        gdal.GDT_Float32)

    dataset.SetGeoTransform(get_geotransform(transformparams))

    dataset.SetProjection(wkt_projection)

    dataset.GetRasterBand(1).WriteArray(data)
    dataset.FlushCache()
    return dataset


def threshold_calc(array, name, threshold, transformparams):
    geotiffData = array.ReadAsArray()

    [rows, cols] = geotiffData.shape

    temp1 = np.less(geotiffData, threshold[0])
    np.putmask(geotiffData, temp1, -9999)

    temp2 = np.greater(geotiffData, threshold[1])
    np.putmask(geotiffData, temp2, -9999)

    temp3 = np.logical_and([np.greater_equal(geotiffData, threshold[0])], [np.less_equal(geotiffData, threshold[1])])
    np.putmask(geotiffData, temp3, 1)

    temp4 = np.equal(geotiffData, -9999)
    np.putmask(geotiffData, temp4, 0)

    outData = driver.Create(
        name,
        cols,
        rows,
        1,
        gdal.GDT_Float32)

    outData.SetGeoTransform(get_geotransform(transformparams))

    outData.SetProjection(wkt_projection)
    outData.GetRasterBand(1).WriteArray(geotiffData)
    outData.GetRasterBand(1).SetNoDataValue(-9999)
    outData = None


def models_threshold(giorno, modello):
    for i in range(3):
        with open(modello + ".gra", 'rb') as f:
            byte = calc_byte(i, 0)
            tmpsfc = read(f, nz, byte)
            threshold_calc(
                write_geotiff_file(
                    tmpsfc,
                    tmp_directory + 'tmpsfc_Run' + str(i) + driver_ext,
                    'model'
                ),
                tmp_directory + 'tmpsfc_Run' + str(i) + '_threshold_' + giorno + driver_ext,
                get_threshold('threshold_tmpsfc'),
                'model'
            )
            dataset = None
            f.close

        with open(modello + ".gra", 'rb') as f:
            byte = calc_byte(i, 2)
            rhsfc = read(f, nz, byte)
            threshold_calc(
                write_geotiff_file(
                    rhsfc,
                    tmp_directory + 'rhsfc_Run' + str(i) + driver_ext,
                    'model'
                ),
                tmp_directory + 'rhsfc_Run' + str(i) + '_threshold_' + giorno + driver_ext,
                get_threshold('threshold_rhsfc'),
                'model'
            )
            dataset = None
            f.close

        with open(modello + ".gra", 'rb') as f:
            byte = calc_byte(i, 3)
            wind10m = read(f, nz, byte)
            wind10m = wind10m * 3.6
            threshold_calc(
                write_geotiff_file(
                    wind10m,
                    tmp_directory + 'wind10m_Run' + str(i) + driver_ext,
                    'model'
                ),
                tmp_directory + 'wind10m_Run' + str(i) + '_threshold_' + giorno + driver_ext,
                get_threshold('threshold_wind10m'),
                'model'
            )
            dataset = None
            f.close


def risk_threshold(giorno, rischio_dir):
    runRange = {'1': 'Run0', '2': 'Run1', '3': 'Run2'}

    for value in runRange.values():
        src_ds_dmc = gdal.Open(rischio_dir + "/" + "modello_Rdmc_" + value + "_" + giorno + ".rst")
        src_band_dmc = src_ds_dmc.GetRasterBand(1).ReadAsArray()
        threshold_calc(
                write_geotiff_file(
                    src_band_dmc,
                    tmp_directory + 'rdmc_' + value + '_' + giorno + driver_ext,
                    'rst'
                ),
            tmp_directory + 'rdmc_' + value + '_threshold_' + giorno + driver_ext,
            get_threshold('threshold_dmc'),
            'rst'
        )
        dataset = None
        src_ds_dmc = None

    for value in runRange.values():
        src_ds_ffmc = gdal.Open(rischio_dir + "/" + "modello_Rfff_" + value + "_" + giorno + ".rst")
        src_band_ffmc = src_ds_ffmc.GetRasterBand(1).ReadAsArray()
        threshold_calc(
            write_geotiff_file(
                src_band_ffmc,
                tmp_directory + 'rfff_' + value + '_' + giorno + driver_ext,
                'rst'
            ),
            tmp_directory + 'rfff_' + value + '_threshold_' + giorno + driver_ext,
            get_threshold('threshold_ffmc'),
            'rst'
        )
        dataset = None
        src_ds_ffmc = None


"""
Calcolo numero giorni senza pioggia
"""
def prec_threshold(giorno):
    def prec_sum(prec_data):
        prec_data_len = len(prec_data)
        prec_num_sum = []
        for i in range(prec_data_len):
            prec = prec_data[i].ReadAsArray()

            temp1 = np.less(prec, 5)
            np.putmask(prec, temp1, 1)

            temp2 = np.greater_equal(prec, 5)
            np.putmask(prec, temp2, -9999)

            temp3 = np.equal(prec, -9999)
            np.putmask(prec, temp3, 0)

            prec_num_sum.append(prec)

        prec_tot = sum(prec_num_sum)

        outData = driver.Create(
            tmp_directory + 'prec_' + giorno + driver_ext,
            239,
            263,
            1,
            gdal.GDT_Float32)

        outData.SetGeoTransform(get_geotransform('prec'))

        outData.SetProjection(wkt_projection_prec)
        outData.GetRasterBand(1).WriteArray(prec_tot)
        outData.GetRasterBand(1).SetNoDataValue(-9999)

        # CREO LE SOGLIE PER LA PIOGGIA
        outDataArray = outData.ReadAsArray()

        temp1 = np.less(outDataArray, get_threshold('threshold_prec')[0])
        np.putmask(outDataArray, temp1, -9999)

        temp2 = np.greater(outDataArray, get_threshold('threshold_prec')[1])
        np.putmask(outDataArray, temp2, -9999)

        temp3 = np.logical_and([np.greater_equal(outDataArray, get_threshold('threshold_prec')[0])],
                               [np.less_equal(outDataArray, get_threshold('threshold_prec')[1])])
        np.putmask(outDataArray, temp3, 1)

        temp4 = np.equal(outDataArray, -9999)
        np.putmask(outDataArray, temp4, 0)

        outDataPrec = driver.Create(
            tmp_directory + 'prec_threshold_' + giorno + driver_ext,
            239,
            263,
            1,
            gdal.GDT_Float32)

        outDataPrec.SetGeoTransform(get_geotransform('prec'))

        outDataPrec.SetProjection(wkt_projection_prec)
        outDataPrec.GetRasterBand(1).WriteArray(outDataArray)
        outDataPrec.GetRasterBand(1).SetNoDataValue(-9999)

        outData = None
        outDataPrec = None

    def daterange(start_date, end_date):
        for n in range(int((end_date - start_date).days)):
            yield start_date + timedelta(n)

    start_date = datetime.strptime(giorno, "%Y-%m-%d").date() # - timedelta(days=1)
    end_date = start_date - timedelta(days=7)
    list_prec = []
    for single_date in daterange(end_date, start_date):
        #print(single_date.strftime("%Y-%m-%d"))
        prec_date = single_date.strftime("%Y-%m-%d")
        src_ds_prec = gdal.Open("../../metout/toscana_Prec_dem1000_1_1_263_239_" + prec_date + ".rst")
        list_prec.append(src_ds_prec)
    src_ds_prec = None
    prec_sum(list_prec)


def tot_threshold(giorno):
    directory = os.fsencode(tmp_directory)
    filelists = os.listdir(directory)
    for i in range(3):
        array_mul = []
        for file in filelists:
            filename = os.fsdecode(file)
            if fnmatch.fnmatch(filename, '*Run' + str(i) + '_threshold_' + giorno + driver_ext) or fnmatch.fnmatch(filename, 'prec_threshold_' + giorno + driver_ext):
                raster = gdal.Open(tmp_directory + filename)
                if fnmatch.fnmatch(filename, 'prec_threshold_' + giorno + driver_ext):
                    ds_array = raster.ReadAsArray()
                    array_mul.append(ds_array)
                else:
                    """Only for gdal 1.9.0"""
                    subprocess.call('gdalwarp ' + tmp_directory + filename + ' ' + tmp_directory + 'temp.rst -s_srs "EPSG:4326" -t_srs "EPSG:32632" -of RST -tr 1000 1000 -r bilinear -te 548174.9000000000232831 4670721.0000000000000000 787174.9000000000232831 4933721.0000000000000000', shell=True)

                    """PROVA"""
                    rasterrst = gdal.Open(tmp_directory + 'temp.rst')
                    ds_array = rasterrst.ReadAsArray()
                    array_mul.append(ds_array)
                    rasterrst = None
                    os.remove(tmp_directory + 'temp.rst')

                    """For gdal > 2"""
                    """BUONO"""
                    #ds = gdal.Warp(
                    #    '',
                    #    raster,
                    #    targetAlignedPixels=False,
                    #    outputBounds=[
                    #        548174.9000000000232831,
                    #        4670721.0000000000000000,
                    #        787174.9000000000232831,
                    #        4933721.0000000000000000
                    #    ],
                    #    dstSRS='EPSG:32632',
                    #    srcSRS='EPSG:4326',
                    #    format='VRT',
                    #    outputType=gdal.GDT_Float32,
                    #    xRes=1000,
                    #    yRes=1000,
                    #    resampleAlg=gdalconst.GRIORA_NearestNeighbour
                    #)
                    #ds_array = ds.ReadAsArray()
                    #array_mul.append(ds_array)
                ds_array = None
                ds = None
                raster = None

        final_data = reduce((lambda x, y: x * y), array_mul)

        array_mul_data = driver.Create(
            '/mnt/hd/operativo/risout_prev/arw/fire_presc_threshold_Run' + str(i) + '_' + giorno + driver_ext,
            239,
            263,
            1,
            gdal.GDT_Byte)

        array_mul_data.SetGeoTransform(get_geotransform('prec'))

        array_mul_data.SetProjection(wkt_projection_prec)
        array_mul_data.GetRasterBand(1).WriteArray(final_data)
        array_mul_data.GetRasterBand(1).SetNoDataValue(0)

        array_mul_data = None

    #for file in filelists:
    #    filename = os.fsdecode(file)
    #    os.remove(tmp_directory + filename)

def print_error_log(day, log, e):
        logger = logging.getLogger('fuoco_prescritto')
        hdlr = logging.FileHandler('/mnt/hd/operativo/log/fuoco_prescritto_{1}_{0}.log'.format(day, log))
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr)
        logger.setLevel(logging.ERROR)
        logger.error('Error: {0}'.format(e))

        sys.exit(2)

def print_usage():
    print("calc_fuoco_prescritto.py -d <day> -m <model> -r <rischio>")


def main(argv):
    #print("ARGV      :", sys.argv[1:])
    try:
        # opts is a list of returning key-value pairs, args is the options left after striped
        # the short options 'hi:o:', if an option requires an input, it should be followed by a ":"
        # the long options 'ifile=' is an option that requires an input, followed by a "="
        opts, args = getopt.getopt(argv, "hd:m:r:", ["help", "day=", "model=", "rischio="])
    except getopt.GetoptError as err:
        print(err)  # will print something like "option -a not recognized"
        print_usage()
        sys.exit(2)
    if not opts:
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("calc_fuoco_prescritto.py -d <day> -m <model> -r <rischio>")
            sys.exit(2)
        elif opt in ("-d", "--day"):
            day = arg
        elif opt in ("-m", "--model"):
            model = arg
        elif opt in ("-r", "--rischio"):
            rischio_dir = arg
        else:
            assert False, "unhandled option"

    try:
        """Funzioni per calcolare le soglie"""
        """
        Calcolo le soglie per il modello ../metprev/incendi_arw_ecm_3km_run00.gra
        Le variabili che mi interessano del modello sono:
         - Temperatura °: threshold_tmpsfc = [0, 18]
         - Umidità %: threshold_rhsfc = [40, 75]
         - Intensità del vento (km hr -1): threshold_wind10m = [1, 15]

        Il modello contiene le variabili con media giornaliera giornaliera
        Contiene 4 tempi -> OGGI, DOMANI, DOPO DOMANI, TERZO GIORNO
        """
        models_threshold(day, model)
    except Exception as e:
        print_error_log(day, 'models_threshold', e)

    try:
        """
        Calcolo le soglie per gli indici di rischio FFMC() e DMC()
        "../risout_prev/modello_Rfff_Run0_2020-01-29.rst"
        "../risout_prev/modello_Rfff_Run1_2020-01-29.rst"
        "../risout_prev/modello_Rfff_Run2_2020-01-29.rst"

        Run0, Run1, Run2
        threshold_ffmc = [80, 95]
        threshold_dmc = [0, 20]
        """
        risk_threshold(day, rischio_dir)
    except Exception as e:
        print_error_log(day, 'risk_threshold', e)

    try:
        """
        Calcolo le soglie per i giorni di pioggia
        "../metout/toscana_Prec_dem1000_1_1_263_239_"+prec_date+".rst"
        threshold_prec = [1, 7]
        """
        prec_threshold(day)
    except Exception as e:
        print_error_log(day, 'prec_threshold', e)

    try:
        """
        Calcolo della maschera finale VERO-FALSO per individuare
        le finestre ambientali per l'applicazione del fuoco prescritto.
        """
        tot_threshold(day)
    except Exception as e:
        print_error_log(day, 'tot_threshold', e)

    # updateAITcache(db_connection, to_data)
    # db_connection.close()

if __name__ == '__main__':
    main(sys.argv[1:])
    print("--- %s seconds ---" % (time.time() - start_time))
