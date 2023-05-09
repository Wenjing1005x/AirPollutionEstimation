# For PM2.5
import arcpy
from arcpy import env
from arcpy.sa import *
from arcpy import analysis
import arcpy.da
import string
import math
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


env.overwriteOutput = True
#env.workspace = 'D:/Project Exercise/testpm/outputs_all_pm10'
# Get parameter from metadata file(img_meta.TXT)
def getmeta(meta_file):
    with open(meta_file) as meta:
        row_list = []
        surface_para = {}
        toa_para = {}
        me = meta.readlines()
        for line in me:
            li_split = line.split(", ")
            for r in li_split:
                r_split = r.split(" = ")
                row_list.append(r_split)
                if r_split[0] == '    SUN_ELEVATION':
                    solar_elevation = float(r_split[1])
                    solar_zenith = (90 - solar_elevation) * math.pi / 180
        # get surface and top of atmosphere rescaling parameters
        for item in row_list:
            if item[0] == '  GROUP' and item[1] == 'LEVEL2_SURFACE_REFLECTANCE_PARAMETERS\n':
                idx_surface_start = row_list.index(item)
            elif item[0] == '  END_GROUP' and item[1] == 'LEVEL2_SURFACE_REFLECTANCE_PARAMETERS\n':
                idx_surface_end = row_list.index(item)
            elif item[0] == '  GROUP' and item[1] == 'LEVEL1_RADIOMETRIC_RESCALING\n':
                idx_toa_start = row_list.index(item)
            elif item[0] == '  END_GROUP' and item[1] == 'LEVEL1_RADIOMETRIC_RESCALING\n':
                idx_toa_end = row_list.index(item)
        for row in row_list[idx_surface_start+1:idx_surface_end]:
            i = 1
            while i <= 7:
                if row[0] == '    REFLECTANCE_MULT_BAND_' + str(i):
                    surface_para['sur_multi_b' + str(i)] = float(row[1])
                elif row[0] == '    REFLECTANCE_ADD_BAND_' + str(i):
                    surface_para['sur_add_b' + str(i)] = float(row[1])
                i += 1
        for row in row_list[idx_toa_start+1:idx_toa_end]:
            i = 1
            while i <= 7:
                if row[0] == '    REFLECTANCE_MULT_BAND_' + str(i):
                    toa_para['toa_multi_b' + str(i)] = float(row[1])
                elif row[0] == '    REFLECTANCE_ADD_BAND_' + str(i):
                    toa_para['toa_add_b' + str(i)] = float(row[1])
                i += 1
    # merge dictionaries of surface and TOA parameters
    meta_dict = {**surface_para, **toa_para}
    meta_dict['sun_zenith'] = solar_zenith
    #meta_dict = {'sur_multi_b1': value1,'sur_multi_b2':value2,...}
    return meta_dict


# Create metadata dictionary and image dictionaries for 3 bands used
metaDict = {}
monitorDict = {}
band2_Dict = {}
band3_Dict = {}
band4_Dict = {}

# input boundary of study area
env.workspace = arcpy.GetParameterAsText(0)
boundary = arcpy.GetParameterAsText(1)
# e.g. Delhi can be mosaic by path147row40 and path146row40, then its NumOfMosaic should be 2.
NumOfMosaic = int(arcpy.GetParameterAsText(2))
OutputFileName = arcpy.GetParameterAsText(3)


# input all images and their metadata files,
# use these path to fill metadata dictionary and 3 image dictionaries
for num in range(1, NumOfMosaic+1):
    metaDict[num] = arcpy.GetParameterAsText(5*num-1)
    monitorDict[num] = arcpy.GetParameterAsText(5*num)
    band2_img = arcpy.GetParameterAsText(5*num+1)
    band3_img = arcpy.GetParameterAsText(5*num+2)
    band4_img = arcpy.GetParameterAsText(5*num+3)
    band2_Dict[num] = band2_img
    band3_Dict[num] = band3_img
    band4_Dict[num] = band4_img
arcpy.AddMessage("Input raster data and metadata completed!")
# Dictionaries should be like this:
# monitorDict = {1:'PM2.5 monitored on date of path146row40', 2:'PM2.5 monitored on date of path147row40', ......}
# band2_Dict = {1: 'band1 image of path146row40', 2:'band1 image of path147row40', 3:'', 4:''....}
# band3_Dict = {1: 'band2 image of path146row40', 2:'band2 image of path147row40', 3:'', 4:''....}
# band4_Dict = {1: 'band5 image of path146row40', 2:'band5 image of path147row40', 3:'', 4:''....}

# Extract all images by boundary shapefile, update the image dictionaries
for key in band2_Dict:
    image = band2_Dict[key]
    ExtractedImg = arcpy.sa.ExtractByMask(image, boundary)
    band2_Dict[key] = ExtractedImg
for key in band3_Dict:
    image = band3_Dict[key]
    ExtractedImg = arcpy.sa.ExtractByMask(image, boundary)
    band3_Dict[key] = ExtractedImg
for key in band4_Dict:
    image = band4_Dict[key]
    ExtractedImg = arcpy.sa.ExtractByMask(image, boundary)
    band4_Dict[key] = ExtractedImg
arcpy.AddMessage("Extract study area completed!")
# After this image dictionaries are updated as:
# band2_Dict = {1: 'band2 study area of path146row40', 2:'band2 study area of path147row40', 3:'', 4:''....}
# band3_Dict = {1: 'band3 study area of path146row40', 2:'band3 study area of path147row40', 3:'', 4:''....}
# band4_Dict = {1: 'band4 study area of path146row40', 2:'band4 study area of path147row40', 3:'', 4:''....}

# Get parameters from metadata files
Ratm_dict = {}
for key in metaDict:
    img_callmeta = metaDict[key]
    meta_para = getmeta(img_callmeta)
    '''me_list = getmeta(img_callmeta)
    radio_para = me_list[0]
    surface_para = me_list[1]
    solar_zenith = me_list[2]'''
    # Calculate Rr: Rr = (multi * DN + add)
    Rr_b2 = (meta_para['sur_multi_b2'] * band2_Dict[key] + meta_para['sur_add_b2'])
    Rr_b3 = (meta_para['sur_multi_b3'] * band3_Dict[key] + meta_para['sur_add_b3'])
    Rr_b4 = (meta_para['sur_multi_b4'] * band4_Dict[key] + meta_para['sur_add_b4'])
    # Calculate Rs: Rs = (multi * DN + add) / cos(solar zenith)
    Rs_b2 = (meta_para['toa_multi_b2'] * band2_Dict[key] + meta_para['toa_add_b2']) / math.cos(meta_para['sun_zenith'])
    Rs_b3 = (meta_para['toa_multi_b3'] * band3_Dict[key] + meta_para['toa_add_b3']) / math.cos(meta_para['sun_zenith'])
    Rs_b4 = (meta_para['toa_multi_b4'] * band4_Dict[key] + meta_para['toa_add_b4']) / math.cos(meta_para['sun_zenith'])
    # Calculate Ratm: Ratm = Rs - Rr
    Ratm_b2 = Rs_b2 - Rr_b2
    Ratm_b3 = Rs_b3 - Rr_b3
    Ratm_b4 = Rs_b4 - Rr_b4
    # Add atmosphere reflectance output to TOA dictionary
    Ratm_dict[key] = [Ratm_b2, Ratm_b3, Ratm_b4]

arcpy.AddMessage("Calculate Ratm completed!")
# arcpy.AddMessage(type(Ratm_dict[1][0]))
# After this, Ratm_dict is: (1-path146row40, 2-path147row40, 3-......)
# Ratm_dict = {1:[Ratm146_b2, Ratm146_b3, Ratm146_b4], 2:[Ratm147_b2, Ratm147_b3, Ratm147_b4], 3:[], 4:[]...}

# monitorDict = {1:'PM2.5 monitored on date of path146row40', 2:'PM2.5 monitored on date of path147row40', ......}
# ExPointDict = {1:[ExtractPt_value_b1_image1,ExtractPt_value_b2_image1,ExtractPt_value_b5_image1],
#                2:[ExtractPt_value_b1_image2,....],...}
ExPointDict = {}
for key in Ratm_dict:
    ExtractValuesToPoints(monitorDict[key], Ratm_dict[key][0], f'ExtractPt_value_b2_{str(key)}')
    ExtractValuesToPoints(monitorDict[key], Ratm_dict[key][1], f'ExtractPt_value_b3_{str(key)}')
    ExtractValuesToPoints(monitorDict[key], Ratm_dict[key][2], f'ExtractPt_value_b4_{str(key)}')
    expoint1= f'ExtractPt_value_b2_{str(key)}'
    expoint2 = f'ExtractPt_value_b3_{str(key)}'
    expoint3 = f'ExtractPt_value_b4_{str(key)}'
    ExPointDict[key] =[expoint1, expoint2, expoint3]
arcpy.AddMessage("Extract raster value completed!")

# select valid point where rastervalue is not null and pm10>0
validPointDict = {}
# the structure of validPointDict is the same as ExPointDict but it has the valid point value in
for key in ExPointDict:
    arcpy.analysis.Select(ExPointDict[key][0], f'Selected_b2_{str(key)}', "RASTERVALU > 0 And PM10 > 0")
    arcpy.analysis.Select(ExPointDict[key][1], f'Selected_b3_{str(key)}', "RASTERVALU > 0 And PM10 > 0")
    arcpy.analysis.Select(ExPointDict[key][2], f'Selected_b4_{str(key)}', "RASTERVALU > 0 And PM10 > 0")
    slcted_point1 = f'Selected_b2_{str(key)}'
    slcted_point2 = f'Selected_b3_{str(key)}'
    slcted_point3 = f'Selected_b4_{str(key)}'
    validPointDict[key] = [slcted_point1, slcted_point2, slcted_point3]


# extract pm10 and rastervalue
def field_to_list(inputshp):
    shpfields = ['PM10', 'RASTERVALU']
    pm10 = []
    ras_value = []
    shprows = arcpy.SearchCursor(inputshp,shpfields)
    while True:
        shprow = shprows.next()
        if not shprow:
            break
        pm10.append(shprow.PM10)
        ras_value.append(shprow.RASTERVALU)
    return pm10, ras_value


reg_dict = {}
# ExPointDict = {1:[ExtractPt_value_b2_image1,ExtractPt_value_b3_image1,ExtractPt_value_b4_image1],
#                2:[ExtractPt_value_b2_image2,....],...}
# reg_dict {1:[[pm10 list of image 1, ras_value band1 list of image 1],
#              [pm10 list of image 1,ras_value band2 list of image 1]],
#          2:[pm10 list of image 2,...],...}
for key in validPointDict:
    # reg_dict[key][0] ia as a tuple
    reg_dict2 = field_to_list(validPointDict[key][0])
    reg_dict3 = field_to_list(validPointDict[key][1])
    reg_dict4 = field_to_list(validPointDict[key][2])
    reg_dict[key] = [reg_dict2, reg_dict3, reg_dict4]


# Define the form of the fit,x_data:independent variable,coef_:the parameters needed to be fitted
def func(x_data, coef_a ,coef_b,coef_c,intercept):
    a = x_data[0][:]
    b = x_data[1][:]
    c = x_data[2][:]
    return coef_a*a+coef_b*b+coef_c*c+intercept
pm10_data = {}
# pm10_data {1:pm2.5 list of image 1,2: pm2.5 list of image 2,.....}
for key in reg_dict:
    # select pm2.5 data which is the same as long as the image is the same regardless of band
    pm10_data[key] = np.array(reg_dict[key][0][0])
band_data = {}
# reg_dict {1:[[pm10 list of image 1, ras_value band1 list of image 1],
#              [pm10 list of image 1,ras_value band2 list of image 1],...],
#           2:[pm10 list of image 2,...],...}
# band_data {1:[raster value list of band 1 image 1,
#               raster value list of band 2 image 1,
#               raster value list of band 5 image 1],
#            2:[raster value list of band 1 image 2,.....],...}
for key in reg_dict:
    band_data2 = np.array(reg_dict[key][0][1])
    band_data3 = np.array(reg_dict[key][1][1])
    band_data4 = np.array(reg_dict[key][2][1])
    band_data[key] = [band_data2, band_data3, band_data4]
# get bands data of one image to be independent variables
data_x = {}
for key in reg_dict:
    data_x[key] = np.array([band_data[key][0],band_data[key][1],band_data[key][2]])
# to fit
coefficients = {}
for key in data_x:
    coefficients[key] = curve_fit(func, data_x[key], pm10_data[key])
# get coef_a ,coef_b,coef_c,intercept value as a list,key represents image
coef_dict = {}
for key in coefficients:
    coef_dict[key] = coefficients[key][0]
# Ratm_dict = {1:[Ratm146_b1, Ratm146_b2, Ratm146_b5],
#              2:[Ratm147_b1, Ratm147_b2, Ratm147_b5], 3:[], 4:[]...}
# Calculate PM10: y = (coef_a*band1)+(coef_b*band2)+(coef_c*band5)+intercept
pm10_cal = {}
for key in coef_dict:
    pm10_cal[key] = coef_dict[key][0]*Ratm_dict[key][0] + coef_dict[key][1]*Ratm_dict[key][1] + coef_dict[key][2]*Ratm_dict[key][2] + coef_dict[key][3]

arcpy.AddMessage("Calculate PM10 completed!")

img_list = []
for key in pm10_cal:
    img_list.append(pm10_cal[key])

addpath = os.path.join(env.workspace, OutputFileName)
os.makedirs(addpath, exist_ok=True)
for i in range(len(img_list)):
    arcpy.RasterToOtherFormat_conversion(img_list[i], addpath, 'TIFF')

# Mosaic PM2.5 output layers
env.workspace = addpath
mos_list = arcpy.ListRasters()
if NumOfMosaic > 1:
    arcpy.Mosaic_management(mos_list[1:], mos_list[0])


