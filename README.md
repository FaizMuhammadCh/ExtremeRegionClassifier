# Extreme-Regions
*Gets the regions and timestamps of regions with extreme conditions using NetCDF files and ERA5 data.*

In extremeRegions.py you will have to specify the path to a directory in your computer that has NetCDF file. Along with the path, you need to give it minimum longitude, maximum longitude, maximum latitude and minimum latitude, in the same order. 

After running the the code, the program will ask you to choose variables that you want to evaluate, for example: Significant wave height, mean wave period etc. This program uses the beaufort scale to get the extreme regions when significant wave height, 10m u-component of wind, or 10m v-component of wind are used. Otherwise, it takes the regions that have the top 5% of highest values.
