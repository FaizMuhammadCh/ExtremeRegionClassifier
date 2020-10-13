!pip install netCDF4
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import numpy as np
from numpy import *
import os
import csv
import sys
%matplotlib inline 

'''
Author: Faiz Muhammad
Description: findExtremeRegions is the main function that calls iterates over all the NetCDF files and extracts extreme regions
by calling other functions. Choose Variable function takes input from the user about which variable to use.
'''

#Input Variables
source = "C:/Users/HP/Documents/Jobs/OSM/Internship/NetCDFfiles/files"
destination = "C:/Users/HP/Documents/Jobs/OSM/Internship/NetCDFfiles/files/output"
waveFileName = "waveResults.txt" #output file 1
windFileName = "windResults.txt"#output file 2
fileCount = 3 #number of latest files you want to select from the source folder
varFlag = 1 #varFlag == 1 when using wave variable and 0 when using wind variable
variable = 'WVDIR_meansealevel' #add the name of your variable
beaufortNum = 6

#defining beaufort scale
bfScale = {}
bfScale[6] = [13.8, 4]
bfScale[7] = [17.1, 5.5]
bfScale[8] = [20.7, 7.5]
bfScale[9] = [24.4, 10]
bfScale[10] = [28.4, 12.5]
bfScale[11] = [32.6, 16]
bfScale[12] = [32.7, 14]
        
def findExtremeRegions(directory, outputLoc, wave, wind, varflag, variable,fileCount,beaufortNum, minLon, maxLon, minLat, maxLat):
        counter = fileCount
        APIinfo = []
        P = []
        bs = 0
        f = open(os.path.join(outputLoc, wind), "a+")
        f.write("Index number; Min Latitude; Max Latitude; Min Longitude; Max Longitude; Time\n")
        
        g = open(os.path.join(outputLoc, wave), "a+")
        g.write("Index number; Min Latitude; Max Latitude; Min Longitude; Max Longitude; Time\n")
        for filename in sorted(os.listdir(directory), reverse = True): #looping through the directory
            if filename.endswith(".nc") and counter > 0:
                print("File Name = ", filename)
                ds, var = chooseVar(os.path.join(directory, filename), variable)
                if ('expver' not in ds.dims): 
                    bs = ds
                    if var != 0:
                        p = ExtremeValuesRegions(ds, varflag, var,beaufortNum, minLon, maxLon, minLat, maxLat) #finding extreme values
                        P.append(p)
                        if (varflag == 1):
                            APIinfo.append(getDetails(p, g, outputLoc)) #appending API information for a dataset
                        elif (varflag == 0):
                            APIinfo.append(getDetails(p, f, outputLoc))
                    else:
                        print("Variable Not Found")
                    counter -= 1
        return APIinfo,P, bs

def chooseVar(path, variable):
        ds = xr.open_dataset(path)
        Flag = False
        for varname in ds:
            if varname == variable:
                Flag = True
                var = varname
                break   
        if Flag == True:
            return ds, var
        else:
            return ds, 0
    
def ExtremeValuesRegions(ds, varflag, varname,bfNum, minLon, maxLon, minLat, maxLat): #calculating extreme values on beaufort 8 or z-vals
        windThresh, waveThresh = getbeaufortScale(bfNum, bfScale)
        var = ds[varname]
        var_denseVessels = var.sel(longitude = slice(minLon, maxLon), latitude = slice (minLat, maxLat)) #need to have global data
        if(varflag == 1):#if the variable selected is significant wave height
            z, threshold = getZvals(var_denseVessels, waveThresh) #getting z values and threshold in z scale for beaufort 8 = 7
        elif(varflag == 0 and varname == "u10"):
            var2 = ds["v10"]
            var2_denseVessels = var2.sel(longitude = slice(minLon, maxLon), latitude = slice (minLat, maxLat))
            magnitudeWind = getMagnitude(var_denseVessels,var2_denseVessels)
            z, threshold = getZvals(magnitudeWind, windThresh) #getting z values and threshold in z scale for beaufort 8 = 20.7
        elif(varflag == 0 and varname == "v10"):
            var2 = ds["u10"]
            var2_denseVessels = var2.sel(longitude = slice(minLon, maxLon), latitude = slice (minLat, maxLat))
            magnitudeWind = getMagnitude(var_denseVessels,var2_denseVessels)
            z, threshold = getZvals(magnitudeWind, windThresh) #getting z values and threshold in z scale for beaufort 8 = 20.7
        else:
            print("Default Condition: No variable found")
            z, trash = getZvals(var_denseVessels, 0)
            threshold = getThreshold_percentile(z, 0.95) #returns the value above which top 5% of z values are
        maskRegion(z, threshold) #converts the regions into binary maps
        return z
    
def getZvals(var_denseVessels, beaufort_threshold):
        z = {}
        #smoothing
        if ('expver' in var_denseVessels.dims):
            var_denseVessels_mean = var_denseVessels.mean(dim = ('time',  'expver')) #computing mean
            var_denseVessels_std = var_denseVessels.std(dim = ('time',  'expver')) #computing standard deviation
            for i in range (0, len(var_denseVessels.time)): #computing z values for each timestamp in a dataset
                if (var_denseVessels[i][0].sum() > 0): 
                    z[i] = (var_denseVessels[i][0] - var_denseVessels_mean)/var_denseVessels_std    
                else:
                    z[i] = (var_denseVessels[i][1] - var_denseVessels_mean)/var_denseVessels_std 
        else:
            var_denseVessels_mean = var_denseVessels.mean(dim = ('time')).values
            var_denseVessels_std = var_denseVessels.std(dim = ('time')).values
            count = 0
            for i in range (0, len(var_denseVessels.time)): #computing z values for each timestamp in a dataset
                if (var_denseVessels[i].sum() > 0): 
                    z[count] = (var_denseVessels[i] - var_denseVessels_mean)/var_denseVessels_std
                    count = count + 1        
        if beaufort_threshold != 0:
            beaufort_threshold = (beaufort_threshold - np.nanmean(var_denseVessels_mean))/np.nanmean(var_denseVessels_std) #beaufort threshold being standardin
        return z, beaufort_threshold
    
def getMagnitude(u10, v10):
        return np.hypot(u10,v10)
    
def getThreshold_percentile(z, top):
        z_total = []
        for i in range (len(z)): 
            z_total.append(z[i].values)
        z_total = np.array(z_total)
        threshold = np.nanquantile(z_total, 0.95)
        return threshold
    
def maskRegion(z, threshold):
        for k in range (len(z)): #len(z)
            result = np.array(z[k])
            where_are_NaNs = isnan(result) 
            result[where_are_NaNs] = 0 #converting NaNs region to 0 as well
            result[result >= threshold] = 100
            result[result < threshold] = 0
            z[k].values = result #storing binary maps in z
            
def getDetails(p, f,outputLoc):
        APIinfo = []
        percentage = 0
        for i in range (0, len(p)):
            if (p[i].sum().values > 0):
                minLat, maxLat, minLon, maxLon = getMinMaxLatLon(p[i])
                if (minLat == -1 and maxLat == -1 and minLon == -1 and maxLon == -1):
                     pass
                else:
                    APIinfo.append((i,minLat,maxLat,minLon,maxLon,p[i].time.values))
                    f.write(str(i)+ "; "+ str(minLat)+ "; "+ str(maxLat)+ "; "+ str(minLon)+ "; "+ str(maxLon)+ "; "+ str(p[i].time.values)+'\n')
            else:
                pass
        return APIinfo

def getMinMaxLatLon(z):
        lat = []
        lon = []
        if hasattr(ds, 'longitude'):
            lenlong = len(z.longitude)
            lenlat = len(z.latitude)
            latsize = z.latitude.size
            lonsize = z.longitude.size
            latvalues = z.latitude.values
            lonvalues = z.longitude.values
            lon0value = z.longitude[0].values
            lon1value = z.longitude[1].values
            lat0value = z.latitude[0].values
            lat1value = z.latitude[1].values
         else:
            lenlong = len(z.lon)
            lenlat = len(z.lat)
            latsize = z.lat.size
            lonsize = z.lon.size
            latvalues = z.lat.values
            lonvalues = z.lon.values
            lon0value = z.lon[0].values
            lon1value = z.lon[1].values
            lat0value = z.lat[0].values
            lat1value = z.lat[1].values
        
        if latsize < 2 or lonsize < 2: #if both latitude and longitude are arrays
            if (latsize < 2 and lonsize < 2 ):
                minLon = lonvalues
                minLat = latvalues
                if(z[0][0] == 100):
                    lat.append(minLat)
                    lon.append(maxLat)
            elif (latsize < 2): # If only longitude is an array
                minLat = latvalues
                minLon = lon0value
                lonDif = abs(minLon - lon1value) # to see the difference between 1st and 2nd index values
                for i in range(0, lenlong):
                    if z[i] == 100:
                        lat.append(minLat)
                        lon.append(minLon + (i*lonDif))
            elif (lonsize < 2): # if only latitude is an array
                
                minLat = lat0value
                minLon = lonvalues
                latDif = abs(minLat - lat1value)
                for i in range(0, lenlong):
                    if z[i] == 100:
                        lat.append(minLat - (i*latDif))
                        lon.append(minLon)      
        else: #if both longitude and latitude are not arrays
            minLon = lon0value
            minLat = lat0value
            latDif = abs(minLat - lat1value)
            lonDif = abs(minLon - lon1value)
            for i in range(0, lenlat):
                for j in range (0, lenlong): 
                    if (z[i][j] == 100):
                        lat.append(minLat - (i*0.5))
                        lon.append(minLon + (j*0.5)) 
        if (lat == [] or lon == []):
            return -1, -1, -1, -1
        else:    
            return min(lat), max(lat), min(lon), max(lon)
        
def getbeaufortScale(number, bfScale):
        if number in range(6, 13):
            return bfScale[number]
        else:
            print("Taking default beaufort value = 8")
            return bfScale[8]
        
#call main function
APIdet, P, bs = findExtremeRegions(source,destination,waveFileName,windFileName,varFlag, variable,fileCount,beaufortNum,-180, 40,-77.5,-20) #minLon, maxLon, maxlat and minLat
