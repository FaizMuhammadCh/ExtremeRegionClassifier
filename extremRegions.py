!pip install netCDF4
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import numpy as np
from numpy import *
import os
%matplotlib inline 
'''
Author: Faiz Muhammad
Description: findExtremeRegions is the main function that calls iterates over all the NetCDF files and extracts extreme regions
by calling other functions. Choose Variable function takes input from the user about which variable to use.
'''
def findExtremeRegions(directory, minLon, maxLon, minLat, maxLat):
        APIinfo = []
        P = []
        for filename in os.listdir(directory): #looping through the directory
            if filename.endswith(".nc"):
                print("File Name = ", filename)
                ds, var = chooseVar(os.path.join(directory, filename))
                bs = ds
                if var != 0:
                    p = ExtremeValuesRegions(ds, var, minLon, maxLon, minLat, maxLat) #finding extreme values
                    P.append(p)
                    APIinfo.append(getDetails(p)) #appending API information for a dataset
        return APIinfo,P, bs

def chooseVar(path):
        while(1):
            ds = xr.open_dataset(path)
            count = 1
            for varname in ds:
                print(count," ",varname)
                count = count + 1
            print(count,"  skip") #giving option to skip the dataset
            value = input("Choose a number:\n")
            print(f'You entered {value}')
            counter = 1
            if (int(value) == count): #selecting the variable
                return ds, 0
            for varname in ds:
                if (counter == int(value)): 
                    var = varname
                    break
                else:
                    counter += 1
            if (int(value) > count or int(value) < 1): #looping if the number entered is invalid
                print("Invalid Number. Please enter a valid number.")
            else:
                break
    
        return ds, var    
    
def ExtremeValuesRegions(ds,varname, minLon, maxLon, minLat, maxLat): #calculating extreme values on beaufort 8 or z-vals
        var = ds[varname]
        var_denseVessels = var.sel(longitude = slice(minLon, maxLon), latitude = slice (minLat, maxLat)) #need to have global data
        if(varname == "swh"):#if the variable selected is significant wave height
            z, threshold = getZvals(var_denseVessels, 7) #getting z values and threshold in z scale for beaufort 8
        elif(varname == "u10"):
            var2 = ds["v10"]
            var2_denseVessels = var2.sel(longitude = slice(minLon, maxLon), latitude = slice (minLat, maxLat))
            magnitudeWind = getMagnitude(var_denseVessels,var2_denseVessels)
            z, threshold = getZvals(magnitudeWind, 20.7) #getting z values and threshold in z scale for beaufort 8
        elif(varname == "v10"):
            var2 = ds["u10"]
            var2_denseVessels = var2.sel(longitude = slice(minLon, maxLon), latitude = slice (minLat, maxLat))
            magnitudeWind = getMagnitude(var_denseVessels,var2_denseVessels)
            z, threshold = getZvals(magnitudeWind, 20.7) #getting z values and threshold in z scale for beaufort 8
        else:
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
                if (var_denseVessels[i][0].sum() > 0): 
                    z[count] = (var_denseVessels[i] - var_denseVessels_mean)/var_denseVessels_std
                    count = count + 1
        if beaufort_threshold != 0:
            beaufort_threshold = (beaufort_threshold - var_denseVessels_mean.mean())/var_denseVessels_std.mean() #beaufort threshold being standardin
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
            
def getDetails(p):
        APIinfo = []
        for i in range (0, len(p)):
            if (p[i].sum().values > 0):
                minLat, maxLat, minLon, maxLon = getMinMaxLatLon(p[i])
                if (minLat == -1 and maxLat == -1 and minLon == -1 and maxLon == -1):
                     print("No Extreme Region")
                else:
                    APIinfo.append((i,minLat,maxLat,minLon,maxLon,p[i].time.values))
                    print ("Index number = ", i, "Min Latitude = ", minLat, "Max Latitude = ", maxLat, "Min Longitude = ", 
                       minLon, "Max Longitude = ", maxLon, "Time = ", p[i].time.values)
            else:
                print("No Extreme Region")
        return APIinfo

def getMinMaxLatLon(z):
        lat = []
        lon = []
        if z.latitude.size < 2 or z.longitude.size < 2: #if both latitude and longitude are arrays
            if (z.latitude.size < 2 and z.longitude.size < 2 ):
                minLon = z.longitude.values
                minLat = z.latitude.values
                if(z[0][0] == 100):
                    lat.append(minLat)
                    lon.append(maxLat)
            elif (z.latitude.size < 2): # If only longitude is an array
                minLat = z.latitude.values
                minLon = z.longitude[0].values
                lonDif = abs(minLon - z.longitude[1].values) # to see the difference between 1st and 2nd index values
                for i in range(0, len(z.longitude)):
                    if z[i] == 100:
                        lat.append(minLat)
                        lon.append(minLon + (i*lonDif))
            elif (z.longitude.size < 2): # if only latitude is an array
                
                minLat = z.latitude[0].values
                minLon = z.longitude.values
                latDif = abs(minLat - z.latitude[1].values)
                for i in range(0, len(z.latitude)):
                    if z[i] == 100:
                        lat.append(minLat - (i*latDif))
                        lon.append(minLon)      
        else: #if both longitude and latitude are not arrays
            minLon = z.longitude[0].values
            minLat = z.latitude[0].values
            latDif = abs(minLat - z.latitude[1].values)
            lonDif = abs(minLon - z.longitude[1].values)
            for i in range(0, len(z.latitude)):
                for j in range (0, len(z.longitude)): 
                    if (z[i][j] == 100):
                        lat.append(minLat - (i*0.5))
                        lon.append(minLon + (j*0.5)) 
        if (lat == [] or lon == []):
            return -1, -1, -1, -1
        else:    
            return min(lat), max(lat), min(lon), max(lon)

APIdet, P, bs = findExtremeRegions("C:/Users/HP/Documents/NetCDFfiles/FOI",200,210,50,40) #minLon, maxLon, maxlat and minLat
