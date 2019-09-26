#%%
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
from scipy.stats import rankdata
from matplotlib.ticker import MultipleLocator, FormatStrFormatter




def time_add(List, tname):
	for i in range(len(List)):
		List[i] = add_times(List[i], tname)
	return List

#read in UKESM ensemble members
fname = '/gws/nopw/j04/aopp/strat_clim_data/CMIP6/UKESM1/historical/zmzw/ZMZW*.nc'
cubes = iris.load(fname)#, callback=callback)
cubes = cubes.extract(iris.Constraint(air_pressure = 10))
UKESM_cubes = time_add(cubes, 'time')

#reanalysis
fname = '/gws/nopw/j04/aopp/strat_clim_data/reanalysis/ERA-5/daily/u_10hpa_daily_era5_1979-2019.nc'
ERA5 = iris.load_cube(fname)#, callback=callback)
ERA5 = add_times(ERA5, 'time')

#ERA reanlysis data
fname = '/gws/nopw/j04/aopp/odimdore/ERA_data/U/U_ERA_10hPa_daily.nc'
U_ERA_cubes = iris.load(fname)
U_ERA = U_ERA_cubes[0]

#load cube of U fields, GC3 and ERA observation data.
def callback(cube, field, filename): 
    # Remove the history attribute. 
    del cube.attributes['history']
    
    
#ERA 40 data
fname = '/gws/nopw/j04/aopp/odimdore/ERA_data/U/U_ERA40/*.nc'
cubes = iris.load(fname, callback=callback)
del cubes[27]
U_ERA40 = cubes.concatenate_cube()
U_ERA40 = add_times(U_ERA40, 'time')

#ERA data from 2015
fname = '/gws/nopw/j04/aopp/odimdore/ERA_data/U/U_ERA_interim2/*.nc'
cubes = iris.load(fname, callback=callback)
U_ERA2 = cubes.concatenate_cube()
U_ERA2 = add_times(U_ERA2, 'time')

#average ERA data over longitudes for ERA data
U_ERA_lon = U_ERA.collapsed('longitude', iris.analysis.MEAN)
U_ERA40_lon = U_ERA40.collapsed('longitude', iris.analysis.MEAN)
U_ERA2_lon = U_ERA2.collapsed('longitude', iris.analysis.MEAN)

#constrain latitude to 60N to find stratospheric polar jet U for ERA data
U_strat_ERA = U_ERA_lon.extract(lat_constraint)
U_strat_ERA2 = U_ERA2_lon.extract(lat_constraint)
U_strat_ERA40 = U_ERA40_lon.extract(lat_constraint)




month = [' ','Jul',' ', 'Aug',' ', 'Sep',' ', 'Oct',' ', 'Nov' ,' ', 'Dec',' ', 'Jan',' ', 'Feb',' ', 'Mar',' ', 'Apr',' ', 'May',' ', 'Jun']

def day_means(U):
	#take day means and stdevs
	U_new = U.collapsed('longitude', iris.analysis.MEAN).extract(iris.Constraint(latitude=60))
	
	daymeans1 = U_new.aggregated_by('day_number', iris.analysis.MEAN)
	daystd1 = U_new.aggregated_by('day_number', iris.analysis.STD_DEV)

	#recentre daily means
	daystd = np.append(daystd1[180:360].data, daystd1[0:180].data)
	daymean = np.append(daymeans1[180:360].data, daymeans1[0:180].data)

	return daymean, daystd

## take day means of UKESM ensemble data
UKESM_mean = np.zeros([9, 360])
UKESM_std = np.zeros([9, 360])
for i in range(9):
	Set = day_means(UKESM_cubes[i])
	UKESM_mean[i,:] = Set[0]
	UKESM_std[i,:] = Set[1]

UK = np.mean(UKESM_mean, axis = 0)
UKstd = np.sqrt(np.sum(UKESM_std**2, axis = 0)/9)


ERA_daymeans1 = U_strat_ERA.aggregated_by('day_number', iris.analysis.MEAN)
U_strat_ERA402 = iris.cube.CubeList([U_strat_ERA40, U_strat_ERA2]).concatenate_cube()
ERA402_daymeans1 = U_strat_ERA402.aggregated_by('day_number', iris.analysis.MEAN)
ERA_daystd1 = U_strat_ERA.aggregated_by('day_number', iris.analysis.STD_DEV)
ERA402_daystd1 = U_strat_ERA402.aggregated_by('day_number', iris.analysis.STD_DEV)
ERA_daymean = np.append(ERA_daymeans1[180:366].data, ERA_daymeans1[0:180].data)
ERA402_daymean = np.append(ERA402_daymeans1[180:366].data, ERA402_daymeans1[0:180].data)
ERA_daymeanc = np.mean([ERA_daymean, ERA402_daymean], axis = 0)
ERA_daystd = np.append(ERA_daystd1[180:366].data, ERA_daystd1[0:180].data)
ERA402_daystd = np.append(ERA402_daystd1[180:366].data, ERA402_daystd1[0:180].data)
ERA_daystdc = np.mean([ERA_daystd, ERA402_daystd], axis = 0)
ERA_daystdc = 0.5*np.sqrt(ERA_daystd**2 + ERA402_daystd**2)
time1 = np.linspace(0,12,len(UK))
time2 = np.linspace(0,12,len(ERA_daymeanc))

fig, ax = plt.subplots(figsize=(8, 5), dpi = 500) 
#plt.plot([0,12], [0,0])
plt.plot(time1, UK, color= 'blue', label = 'UKESM')
plt.plot(time2, ERA_daymeanc, color= 'black', label = 'ERA-interim')
plt.plot(time1, UK + UKstd, linestyle = '--', color = 'blue')
plt.plot(time1, UK - UKstd, linestyle = '--', color = 'blue')
plt.plot(time2, ERA_daymeanc + ERA_daystdc, linestyle = '--', color = 'black')
plt.plot(time2, ERA_daymeanc - ERA_daystdc, linestyle = '--', color = 'black')

plt.xticks(np.arange(0,12,0.5))
minorLocator = MultipleLocator(1)
ax.set_xticklabels(month)
ax.xaxis.set_minor_locator(minorLocator)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='major',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True)
plt.ylabel('ZMZW (ms$^{-1}$)')
plt.legend()
plt.xlim(0,12)
#plt.title('ZMZW daily climatology, 60N 10hPa')
plt.show()
fig.savefig('../figures/60N_10hPa_for_ACSIS.png', dpi = 200)

