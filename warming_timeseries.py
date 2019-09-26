import iris
import iris.coord_categorisation as coord_cat
import iris.plot as iplt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from GC3_diagnostics_functions import *
import math
import pandas as pd
from collections import Counter
from netCDF4 import num2date, date2num
from datetime import datetime, timedelta
import iris_grib
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

################################################READING IN DATA AND INITIALISING ARRAYS
#start and end years for data
start_yearN96 = 1850
end_yearN96 = 2349
years_N96 = np.arange(start_yearN96, end_yearN96, 1)

#start and end years for data
start_yearUKESM = 2060
end_yearUKESM = 2459
years_UKESM = np.arange(start_yearUKESM, end_yearUKESM, 1)

start_yearN216 = 1850
end_yearN216 = 2149
years_N216 = np.arange(start_yearN216, end_yearN216, 1)

start_yearERA = 1979 
end_yearERA = 2014
years_ERA = np.arange(start_yearERA, end_yearERA, 1)

start_yearERA40 = 1959 
end_yearERA40 = 1978
years_ERA = np.arange(start_yearERA40, end_yearERA40, 1)

alt_constraint = iris.Constraint(Pressure = [10])#,30,50,70,100])

#load cube of U fields, GC3 and ERA observation data.
def callback(cube, field, filename): 
    # Remove the history attribute. 
    del cube.attributes['history']
    del cube.attributes['valid_max']
    del cube.attributes['valid_min']

#GC3 U fields
fname = '../long_runs/CMIP_data/U/N96/daily/u-ar766_U_zonmean10hPa_1850-2349.nc'
cubes = iris.load(fname, lat_constraint)#, callback=callback)
U_N96 = cubes[0]
U_N96 = add_times(U_N96, 'time')
U_N96 = U_N96[:,0,0]

fname = '../long_runs/CMIP_data/U/N216/daily/u-aq281_U10hPazonmean_1850-2349.nc'
cubes = iris.load(fname, lat_constraint)#, callback = callback)
#cubes[-1] = cubes[-1].extract(iris.Constraint(Pressure = 10.))
U_N216 = cubes[0]
U_N216 = add_times(U_N216, 'time')
U_N216 = U_N216[:,0,0]

#UKESM U field
fname = '../long_runs/CMIP_data/UKESM1/U/daily/*.nc'
cubes = iris.load(fname, lat_constraint, callback=callback)
U_UKESM = cubes.concatenate_cube()
U_UKESM = add_times(U_UKESM, 't')
U_UKESM = U_UKESM[:,0,0]

def area_series(U_strat, start_year, end_year):
	U_SSW_seasons, SSW, final, indexi, indexj, SSW_year = [], [], [], [], [], []
	areas = np.empty(0)
	U_strat_ND = U_strat.extract(iris.Constraint(month = ['Nov', 'Dec']))
	U_strat_JFMA = U_strat.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr']))
	for i in range(start_year, end_year-1):
		print i
		dummy2 = iris.cube.CubeList([U_strat_ND.extract(iris.Constraint(year = i)),U_strat_JFMA.extract(iris.Constraint(year = i+1))])
		U_SSW_seasons.append(dummy2.concatenate_cube())
		print U_SSW_seasons
	for j in range(len(U_SSW_seasons)):
		plt.plot(np.arange(180), U_SSW_seasons[j].data)
		plt.show()
		print j	
		i=0
		while i < len(U_SSW_seasons[j].data):
			if U_SSW_seasons[j].data[i] < 0 and U_SSW_seasons[j].data[i-1] >= 0:
				final = i
				print 'here3'
				for x in range(i, len(U_SSW_seasons[j].data)-10):
					if np.all(U_SSW_seasons[j].data[x:x+10] > 0):
						print 'here2'
						i = i+19
					elif x == len(U_SSW_seasons[j].data)-10:
						print 'here3'
						i = 10000000	
					else:
						j=j	
			else:
				i = i+1
		winter = U_SSW_seasons[j].data[0:final]
		area = sum(winter[winter <= 0])*-1				
		areas = np.append(areas, area)
	return areas

areas_N96 = area_series(U_N96, 1851, 1860)


	
























