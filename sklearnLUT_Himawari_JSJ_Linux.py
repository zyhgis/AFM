# -*- coding: utf-8 -*-
import gdal
import os
import time
import math
import numpy as np
from osgeo import gdal
import string
from sklearn.neighbors import KDTree
import pandas as pd
import multiprocessing as mp
from scipy import spatial as sp
import datetime
import sys
class Tiff:
    def read_tif(self, filename):
        dataset = gdal.Open(filename)
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        im_bands = dataset.RasterCount

        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)

        del dataset

        return im_proj, im_geotrans, im_data

    def write_tif(self, filename, im_proj, im_geotrans, im_data):
        if 'int8' in im_data.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in im_data.dtype.name:
            datatype = gdal.GDT_Int16
        else:
            datatype = gdal.GDT_Float32

        if len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape
        else:
            im_bands, (im_height, im_width) = 1, im_data.shape

        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)

        dataset.SetGeoTransform(im_geotrans)
        dataset.SetProjection(im_proj)

        if im_bands == 1:
            dataset.GetRasterBand(1).WriteArray(im_data)
        else:
            for i in range(im_bands):
                dataset.GetRasterBand(i + 1).WriteArray(im_data[i])

        del dataset


def listdir(path):
    list_name = []
    for file in os.listdir(path):
        # print file
        if os.path.splitext(file)[1] == '.tif':
            file_path = os.path.join(path, file)
            list_name.append(file_path)
    return list_name


def ReadLUT(pathLUT):

    lutdic={}
    linelist=open(pathLUT)
    for line in open(pathLUT):
        templut = []
        # rapp=line
        p = ','.join(line.split())
        temp = p.split(",")
        keyvalue = temp[0] + temp[1]
        tempvalue=temp[2:]
        tempvalue=[float(f) for f in tempvalue]
        templut.append(tempvalue)
        if keyvalue in lutdic.keys():
            lutdic[keyvalue].append(tempvalue)
        else:
            lutdic[keyvalue] = templut
    return(lutdic)



def SearchNearLAI(LAI, LAIvaluelist):
    new = [abs(i - LAI) for i in LAIvaluelist]
    a = new.index(min(new))
    return LAIvaluelist[a]


def Search2DLUT(b3valueinput, b4value, AA):
    length = len(AA)
    #print(str(length))
    # AA = lut
    indiceslist = []
    for b3value in [b3valueinput[1]]:
        query = [[b3value, b4value]]
        querylut = np.zeros(shape=(length, 2))
        querylut[:, 0] = AA[:, 20] / query[0][0]
        querylut[:, 1] = AA[:, 24] / query[0][1]
        tree = sp.cKDTree(np.array(querylut))
        dists, indices = tree.query([[1, 1]], k=100)
        indiceslist = indiceslist + list(indices[0])
    # dists, indices = tree.query(query, k=100)
    # print dists
    templailist = []
    tempFPAR1list = []
    tempFPAR2list = []
    tempFPAR3list = []
    # print(len(indiceslist))

    for ii in indiceslist:
        if AA[ii][9]>0.1:
            templailist.append(AA[ii][2])
            tempFPAR1list.append(AA[ii][8])
            tempFPAR2list.append(AA[ii][9])
            tempFPAR3list.append(AA[ii][10])
        # tempFPAR4list.append(templut[ii][14])
        # templailist.append(templut[ii][0])
    datalen=len(templailist)
    if datalen<2:
        tempFAPAR1=0
        templai=0
        tempFAPAR2=0
        tempFAPAR3=0
        laistd=0
        FAPAR1std=0
        FAPAR2std=0
        FAPAR3std=0
    else:
        #templailist.sort()
        #tempFPAR1list.sort()
        #tempFPAR2list.sort()
        #tempFPAR3list.sort()
        templai = np.mean(templailist)
        laistd = np.std(templailist)
        tempFAPAR1 = np.mean(tempFPAR1list)  # difuse
        FAPAR1std = np.std(tempFPAR1list)
        tempFAPAR2 = np.mean(tempFPAR2list)  # direct
        FAPAR2std = np.std(tempFPAR2list)
        tempFAPAR3 = np.mean(tempFPAR3list)  # total
        FAPAR3std = np.std(tempFPAR3list)

    # resultLAI[i][j] = templai*100
    # resultFAPAR[i][j] = tempFAPAR * 100
    return templai, tempFAPAR1, tempFAPAR2, tempFAPAR3, laistd, FAPAR1std, FAPAR2std, FAPAR3std


def RetrivalFAPAR(MDHMlist):
    step = 0.2
    s = os.sep
    basePATH = "/home/ubuntu/Data"
    lutbasePATH = "/home/ubuntu"
    basePATH = "/home/zyh/FAPAR/InputData/"
    lutbasePATH = "/home/zyh/FAPAR/InputData/"
    H8AODBasePATH = basePATH + s + "H8_AOT"
    MODISAODBasePATH = basePATH + s + "MODIS_AOT"

    Band3LUT = basePATH + s + "Band3LUTSKYL.txt"
    Band4LUT = basePATH + s + "Band4LUTSKYL.txt"
    Band3LUTdict = ReadLUT(Band3LUT)
    Band4LUTdict = ReadLUT(Band4LUT)
    bandbasepath = basePATH + s + "Reflectance"
    outputpath = basePATH + s + "FAPAR"
    LAIpath = basePATH + s + "LAI"
    A = Tiff() 
    # MDHMlist = []
    ALLLAIdatadic = {}
    LUTSOZ = [25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70]
    LAIvaluelist = [0.01, 1.5, 2.5, 3.25, 3.75, 4, 4.25, 4.75, 5.5, 6.5]
    lutsoz = 0
    skyllist=np.arange(0.01,0.6,0.1)

    for MDHMHM in MDHMlist:
        MDHM=MDHMHM[0:4]
        hour=MDHMHM[4:6]
        minute=MDHMHM[6:7]

        try:
            H8AOTfile = H8AODBasePATH + s + "JSJ_H08_2019" + MDHM + "_" + hour + "00_1HARP030_FLDK.02401_02401rged.tif"
            if os.path.exists(H8AOTfile)==False:
                continue
            AODdata1 = A.read_tif(H8AOTfile)[2]

            dateday = datetime.datetime(2019, int(MDHM[0:2]), int(MDHM[2:]))
            doy = str(dateday.timetuple().tm_yday).zfill(3)
            modisaodfilelist = listdir(MODISAODBasePATH)
            modisaodlist = []
            for i in modisaodfilelist:
                # print(os.path.split(i)[1][0:30])
                if os.path.split(i)[1][0:30] == "JSJ_Re_ProjectMOD04_L2A2019" + doy or os.path.split(i)[1][0:30] == "JSJ_Re_ProjectMYD04_L2A2019" + doy:
                    modisaodlist.append(i)

            Angleinputpath = basePATH + s + "Himawari8" + s + "2019" + MDHM[0:2] + s + MDHM[2:] + s + "JSJ"

            bandinputpath = bandbasepath + s + "2019" + MDHM[0:2] + s + MDHM[2:] + s + "LUTResult"

            band3name = bandinputpath + s + "Band3" + s + "LbertJSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_03.tif"
            if os.path.exists(band3name) == False:
                continue
            LAIname = LAIpath + s + "JSJNC_H08_2019" + MDHM + "_00_R21_FLDK.06001_06001_LAI.tif"
            if os.path.exists(LAIname) == False:
                continue
            laidata = A.read_tif(LAIname)[2]

            band4name = bandinputpath + s + "Band4" + s + "LbertJSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_04.tif"
            print(band3name)
            # JSJNC_H08_20190801_0010_R21_FLDK.06001_06001_AZ.tif
            sozname = Angleinputpath + s + "SOZ" + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_OZ.tif"
            soaname = Angleinputpath + s + "SOA" + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_OA.tif"

            b3 = A.read_tif(band3name)
            b4 = A.read_tif(band4name)

            b3valuelist = b3[2]
            b4valuelist = b4[2]
            datasize = b3valuelist.shape
            modisAODdata = []
            if len(modisaodlist) > 0:
                for eachmodis in modisaodlist:
                    MODISAOTfile = eachmodis
                    modisAODdata.append(A.read_tif(MODISAOTfile)[2])
            else:
                modisAODdata = [np.ones(datasize) * (-1)]


            soz = A.read_tif(sozname)
            soa = A.read_tif(soaname)
            sozvaluelist = soz[2]
            soavaluelist = soa[2]

            tuple01 = b3valuelist.shape
            # tuple03 = b3value.shape

            row01 = tuple01[0]  
            col01 = tuple01[1]
            sozvaluepure = sozvaluelist[sozvaluelist > 0]
            meansozvalue = np.mean(sozvaluepure) / 100

            yushu = meansozvalue % 10
            if yushu < 2.5:
                SOZ = int(meansozvalue / 10) * 10
            elif yushu < 7.5 and yushu >= 2.5:
                SOZ = int(meansozvalue / 10) * 10 + 5
            elif yushu >= 7.5:
                SOZ = int(meansozvalue / 10) * 10 + 10

            soavaluepure = soavaluelist[soavaluelist > -18100]
            soavaluepure2 = soavaluepure[soavaluepure < 18100]
            meansoavalue = np.mean(soavaluepure2) / 100
            if meansoavalue < 0:
                meansoavalue = 360 + meansoavalue

            if int(meansoavalue / 10) == round(meansoavalue / 10.0, 0):
                SOA = str(int(meansoavalue / 10) * 10)
            else:
                SOA = str(int(meansoavalue / 10) * 10 + 10)
            if SOZ>80:
                SOZ=80
            band4templist = np.array(Band4LUTdict[str(SOZ) + str(70)])
            band3templist = np.array(Band3LUTdict[str(SOZ) + str(70)])

            tempsoz = str(int(SearchNearLAI(meansozvalue, LUTSOZ) * 10))
            if tempsoz != lutsoz:
                lutpath = lutbasePATH + s + "HongheLUT/HongheLUTSOZ" + tempsoz + ".txt"
                lutdata = pd.read_table(lutpath, sep=' ', header=None)
                lutsoz = tempsoz

            resultFAPAR1 = np.zeros((row01, col01), dtype=np.int16)
            resultFAPAR2 = np.zeros((row01, col01), dtype=np.int16)
            resultFAPAR3 = np.zeros((row01, col01), dtype=np.int16)

            # resultLAIstd = np.zeros((row01, col01), dtype=np.int16)
            resultFAPAR1std = np.zeros((row01, col01), dtype=np.int16)
            resultFAPAR2std = np.zeros((row01, col01), dtype=np.int16)
            resultFAPAR3std = np.zeros((row01, col01), dtype=np.int16)

            for i in range(row01):
                # print i
                for j in range(col01):
                    # b1value, b2value, b3value, b4value = 0.1978, 0.1617, 0.1015, 0.5050
                    b3value = b3valuelist[i][j] / 10000.0
                    if b3value <= 0:
                        continue

                    H8AOT = AODdata1[i][j] * 0.0002
                    try:
                        MODISADTlist = []
                        for each in modisAODdata:
                            if each[i][j] > 0:
                                MODISADTlist.append(each[i][j] * 0.0010000000474974513)
                        if len(MODISADTlist)>0:
                            MODISADT = np.mean(MODISADTlist)
                        else:
                            MODISADT=-1
                    except:
                        MODISADT = -1

                    if H8AOT > 0 and H8AOT < 5:
                        AOT = H8AOT
                    elif MODISADT > 0:
                        AOT = MODISADT
                    else:
                        continue
                    AOTlist = band4templist[:, 0]
                    AOTlist2 = np.abs(AOTlist - AOT).tolist()
                    minindex = AOTlist2.index(min(AOTlist2))
                    skyl1 = 1 - band4templist[minindex][1]

                    AOTlist = band3templist[:, 0]
                    AOTlist2 = np.abs(AOTlist - AOT).tolist()
                    minindex = AOTlist2.index(min(AOTlist2))
                    skyl2 = 1 - band4templist[minindex][1]

                    skyl = (skyl1 + skyl2) / 2

                    b3valueinput = [b3value * 0.9, b3value, b3value * 1.1]
                    b4value = b4valuelist[i][j] / 10000.0
                    soavalue = soavaluelist[i][j] / 100.0
                    templai = laidata[i][j] / 1000.0

                    if soavalue < 0:
                        soavalue = 360 + soavalue
                    psi = abs(soavalue - 170)
                    if psi > 180:
                        psi = psi - 180
                    psi = int(psi / 10) * 10
                    data1 = lutdata.loc[(lutdata[1] == psi)]
                    templai = SearchNearLAI(templai, LAIvaluelist)

                    # LAIlist = data1.loc[(data1[2] > templai * (1 - step))]
                    # LAIlist2 = LAIlist.loc[(LAIlist[2] < templai * (1 + step))]
                    LAIlist = data1.loc[(data1[2] > templai-1.2)]
                    LAIlist2 = LAIlist.loc[(LAIlist[2] < templai+1.2)]

                    band4list = LAIlist2.loc[(LAIlist2[24] > b4value * (1 - step))]
                    band4list2 = band4list.loc[(band4list[24] < b4value * (1 + step))]

                    listskyl = SearchNearLAI(skyl, skyllist)
                    # listskyl = SearchNearLAI(skyl, skyllist)
                    band4list = band4list2.loc[(band4list2[7] > listskyl - 0.21)]
                    band4list2 = band4list.loc[(band4list[7] < listskyl + 0.21)]
                    # inputlut = np.array(band4list2.values)
                    inputlut = np.array(band4list2.values)

                    if len(inputlut) <= 100:
                        # LAIlist=band4list2.loc[(band4list2[2] > templai * (1-step-0.1))]
                        # LAIlist2=LAIlist.loc[(LAIlist[2]<templai*(1+step+0.1))]
                        inputlut = np.array(LAIlist2.values)
                        if len(inputlut) <= 100:
                            inputlut = np.array(data1.values)
                            # print(b3valueinput)
                    result = Search2DLUT(b3valueinput, b4value, inputlut)
                    resultFAPAR1[i][j] = result[1] * 1000
                    resultFAPAR2[i][j] = result[2] * 1000
                    # resultFAPAR3[i][j] = result[3] * 1000
                    resultFAPAR3[i][j] = (result[1]*skyl+result[2]*(1-skyl))*1000

                    resultFAPAR1std[i][j] = result[5] * 1000
                    resultFAPAR2std[i][j] = result[6] * 1000
                    resultFAPAR3std[i][j] = result[7] * 1000
            A.write_tif(
                outputpath + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_FAPAR1.tif",
                b3[0], b3[1], resultFAPAR1)  
            A.write_tif(
                outputpath + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_FAPAR2.tif",
                b3[0], b3[1], resultFAPAR2)  
            A.write_tif(
                outputpath + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_FAPAR3.tif",
                b3[0], b3[1], resultFAPAR3) 
            
        except:
            print("Unexpected error:", sys.exc_info()[0]+band3name)


def RetrivalLAI(MDHMlist):
    step = 0.2
    s = os.sep
    # bandinputpath = basePATH+s+"Reflectance"
    basePATH = "/home/ubuntu/Data"
    lutbasePATH = "/home/ubuntu"
    basePATH = "/home/zyh/FAPAR/InputData/"
    lutbasePATH = "/home/zyh/FAPAR/InputData/"
    # basePATH=r"C:\Users\zyh\Documents\FAPAR"

    bandbasepath = basePATH + s + "Reflectance"
    outputpath = basePATH + s + "LAI"
    A = Tiff()  
    # MDHMlist = []
    ALLLAIdatadic = {}
    LUTSOZ = [25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70]
    LAIvaluelist = [0.01, 1.5, 2.5, 3.25, 3.75, 4, 4.25, 4.75, 5.5, 6.5]
    lutsoz = 0
    for MDHM in MDHMlist:
        Angleinputpath = basePATH + s + "Himawari8" + s + "2019" + MDHM[0:2] + s + MDHM[2:] + s + "JSJ"
        for hour in range(24):
            hour = str(hour).zfill(2)
            for minute in range(6):
                minute = str(minute)
                try:
                    bandinputpath = bandbasepath + s + "2019" + MDHM[0:2] + s + MDHM[2:] + s + "LUTResult"

                    band3name = bandinputpath + s + "Band3" + s + "LbertJSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_03.tif"
                    if os.path.exists(band3name) == False:
                        continue

                    band4name = bandinputpath + s + "Band4" + s + "LbertJSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_04.tif"
                    print(band3name)
                    # JSJNC_H08_20190801_0010_R21_FLDK.06001_06001_AZ.tif
                    sozname = Angleinputpath + s + "SOZ" + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_OZ.tif"
                    soaname = Angleinputpath + s + "SOA" + s + "JSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_OA.tif"

                    b3 = A.read_tif(band3name)
                    b4 = A.read_tif(band4name)

                    b3valuelist = b3[2]
                    b4valuelist = b4[2]
                    datasize = b3valuelist.shape

                    ALLLAIdatadic[hour + minute] = np.zeros(datasize, dtype=float, order='C')

                    soz = A.read_tif(sozname)
                    soa = A.read_tif(soaname)
                    sozvaluelist = soz[2]
                    soavaluelist = soa[2]

                    tuple01 = b3valuelist.shape

                    row01 = tuple01[0]  
                    col01 = tuple01[1]
                    sozvaluepure = sozvaluelist[sozvaluelist > 0]
                    meansozvalue = np.mean(sozvaluepure) / 100
                    tempsoz = str(int(SearchNearLAI(meansozvalue, LUTSOZ) * 10))
                    if tempsoz != lutsoz:
                        lutpath = lutbasePATH + s + "HongheLUT/HongheLUTSOZ" + tempsoz + ".txt"
                        lutdata = pd.read_table(lutpath, sep=' ', header=None)
                        # lutdata = ""
                        lutsoz = tempsoz
                    resultLAI = np.zeros((row01, col01), dtype=np.int16)
                    resultLAIstd = np.zeros((row01, col01), dtype=np.int16)

                    for i in range(row01):
                        # print i
                        for j in range(col01):
                            # b1value, b2value, b3value, b4value = 0.1978, 0.1617, 0.1015, 0.5050
                            b3value = b3valuelist[i][j] / 10000.0
                            if b3value <= 0:
                                continue
                            b3valueinput = [b3value * 0.9, b3value, b3value * 1.1]
                            b4value = b4valuelist[i][j] / 10000.0
                            soavalue = soavaluelist[i][j] / 100

                            if soavalue < 0:
                                soavalue = 360 + soavalue
                            psi = abs(soavalue - 170)
                            if psi > 180:
                                psi = psi - 180
                            psi = int(psi / 10) * 10
                            data1 = lutdata.loc[(lutdata[1] == psi)]

                            band4list = data1.loc[(data1[24] > b4value * (1 - step))]
                            band4list2 = band4list.loc[(band4list[24] < b4value * (1 + step))]
                            inputlut = np.array(band4list2.values)

                            if len(inputlut) < 100:
                                inputlut = np.array(data1)
                            # print(str(len(inputlut)))
                            result = Search2DLUT(b3valueinput, b4value, inputlut)
                            # resultLAI[i][j] = result[0] * 1000
                            # resultLAIstd[i][j] = result[4] * 1000
                            ALLLAIdatadic[hour + minute][i][j] = result[0] * 1000
                except:
                    print("Unexpected error:", sys.exc_info()[0]+band3name)


        PixelLAIdatadic = np.zeros(datasize, dtype=float, order='C')
        PixelLAIdatadicstd = np.zeros(datasize, dtype=float, order='C')

        for cloumn in range(datasize[0]):
            for row in range(datasize[1]):
                pixellailist = []
                for key in ALLLAIdatadic.keys():
                    if ALLLAIdatadic[key][cloumn][row] > 0:
                        pixellailist.append(ALLLAIdatadic[key][cloumn][row])
                if len(pixellailist) > 3:
                    pixelmeanlai = np.mean(pixellailist)
                    pixelstdlai = np.std(pixellailist)
                else:
                    pixelmeanlai = 0
                    pixelstdlai = -1
                PixelLAIdatadic[cloumn][row] = pixelmeanlai
                PixelLAIdatadicstd[cloumn][row] = pixelstdlai

        A.write_tif(outputpath + s + "JSJNC_H08_2019" + MDHM + "_00_R21_FLDK.06001_06001_LAI.tif", b3[0], b3[1],
                    PixelLAIdatadic)  
        A.write_tif(outputpath + s + "JSJNC_H08_2019" + MDHM + "_00_R21_FLDK.06001_06001_LAIstd.tif", b3[0], b3[1],
                    PixelLAIdatadicstd)


if __name__ == "__main__":
    s = os.sep
    step = 0.2
    start_time = time.time()
    print(start_time)
    basePATH = "/home/ubuntu/Data"
    lutbasePATH = "/home/ubuntu"
    basePATH = "/home/zyh/FAPAR/InputData/"
    lutbasePATH = "/home/zyh/FAPAR/InputData/"
    bandbasepath = basePATH + s + "Reflectance"

    MDHMlist = []

    for month in range(6, 9):
        MONTH = str(month).zfill(2)
        for day in range(31):
            DAY = str(day).zfill(2)
            MDHM = MONTH + DAY
            nctemp = bandbasepath + s + "2019" + MONTH + s + DAY
            if os.path.exists(nctemp):
                LAIfile=basePATH+"LAI"+s+"JSJNC_H08_2019"+MDHM+"_00_R21_FLDK.06001_06001_LAI.tif"
                print(LAIfile)
                if os.path.exists(LAIfile)==False:
                    print(MDHM)
                    MDHMlist.append(MDHM)
        
    datalen = len(MDHMlist)

    core_num = 16
    start_t = time.time()
    p = mp.Pool(processes=core_num)
    for ii in range(core_num):
       stacksize = int(len(MDHMlist) / core_num)
       tem_file_name = MDHMlist[ii::core_num]
       p.apply_async(RetrivalLAI, args=(tem_file_name,))
    p.close()
    p.join()

    MDHMlist = []
    for hour in range(24):
       hour = str(hour).zfill(2)
       for minute in range(6):
           minute = str(minute)
           for month in range(6, 9):
               MONTH = str(month).zfill(2)
               for day in range(31):
                   DAY = str(day).zfill(2)
                   MDHM = MONTH + DAY
                   MDHMHM = MONTH + DAY + hour + minute
                   bandinputpath = bandbasepath + s + "2019" + MDHM[0:2] + s + MDHM[2:] + s + "LUTResult"
                   band3name = bandinputpath + s + "Band3" + s + "LbertJSJNC_H08_2019" + MDHM + "_" + hour + minute + "0_R21_FLDK.06001_06001_03.tif"
                   FAPARname=basePATH+"FAPAR"+s+"JSJNC_H08_2019"+MDHM+"_"+hour+minute+"0_R21_FLDK.06001_06001_FAPAR1.tif"
                   # JSJNC_H08_20190620_0000_R21_FLDK.06001_06001_FAPAR1.tif
                   if os.path.exists(band3name):
                       if os.path.exists(FAPARname)==False:
                            print("ZhangYinghui"+FAPARname)
                            MDHMlist.append(MDHMHM)
    
    p = mp.Pool(processes=core_num)
    for ii in range(core_num):
       stacksize = int(len(MDHMlist) / core_num)
       tem_file_name = MDHMlist[ii::core_num]
       p.apply_async(RetrivalFAPAR, args=(tem_file_name,))
    p.close()
    p.join()

    print("Search4DLUTrunning time: %s" % ((time.time() - start_time) / 60))
    print("end")
