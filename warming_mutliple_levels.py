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
end_yearN96 = 2351
years_N96 = np.arange(start_yearN96, end_yearN96, 1)
alt_constraint = iris.Constraint(air_pressure = [10,30,50,70,100])

#GC3 U fields
fname = '../long_runs/CMIP_data/U/N96/daily/u-ar766_U_z*.nc'
cubes = iris.load(fname, lat_constraint)
U_N96 = cubes[1][:,:,0]
U_N96 = add_times(U_N96, 'time')
U_N96 = U_N96.extract(alt_constraint)


#SSW_N96, SSW_year_N96, indexi96, indexj96 = SSW_counter(U_N96_10, years_N96[0], years_N96[-1], 0)



con1 = iris.Constraint(month = ['Nov', 'Dec'])
con2 = iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr'])
p = 1
for i in years_N96:
	U = []
	print i
	dummy2 = iris.cube.CubeList([U_N96.extract(iris.Constraint(year = i) & con1) ,U_N96.extract(iris.Constraint(year = i+1) & con2)])
	U.append(dummy2.concatenate_cube())
	print U
	fig = plt.figure(figsize=(18, 12))
	fig.suptitle('ZMZW 10hPa 60N' + str(i) + '/' + str(i+1))
	gs = gridspec.GridSpec(5, 1)
	plt1 = plt.subplot(gs[0,0])
	plt.plot(np.arange(180), U[0][:,4].data)
	plt.plot([0,180], [0,0])
	plt.ylim(-25, 50)
	plt.xlim(0, 180)
	plt.xticks([0,30,60,90,120,150],['01/nov', '01/dec', '01/jan', '01/feb', '01/mar', '01/apr'])
	plt2 = plt.subplot(gs[1,0])
	plt.plot(np.arange(180), U[0][:,3].data)
	plt.plot([0,180], [0,0])
	plt.title('30 hPa')
	plt.ylim(-25, 50)
	plt.xlim(0, 180)
	plt.xticks([0,30,60,90,120,150],['01/nov', '01/dec', '01/jan', '01/feb', '01/mar', '01/apr'])
	plt3 = plt.subplot(gs[2,0])
	plt.plot(np.arange(180), U[0][:,2].data)
	plt.plot([0,180], [0,0])
	plt.title('50 hPa')
	plt.ylim(-25, 50)
	plt.xlim(0, 180)
	plt.xticks([0,30,60,90,120,150],['01/nov', '01/dec', '01/jan', '01/feb', '01/mar', '01/apr'])
	plt4 = plt.subplot(gs[3,0])
	plt.plot(np.arange(180), U[0][:,1].data)
	plt.plot([0,180], [0,0])
	plt.title('70 hPa')
	plt.ylim(-25, 50)
	plt.xlim(0, 180)
	plt.xticks([0,30,60,90,120,150],['01/nov', '01/dec', '01/jan', '01/feb', '01/mar', '01/apr'])
	plt5 = plt.subplot(gs[4,0])
	plt.plot(np.arange(180), U[0][:,0].data)
	plt.plot([0,180], [0,0])
	plt.title('100 hPa')
	plt.ylim(-25, 50)
	plt.xlim(0, 180)
	plt.xticks([0,30,60,90,120,150],['01/nov', '01/dec', '01/jan', '01/feb', '01/mar', '01/apr'])
	plt.tight_layout()
	fig.savefig('../figures/warming_on_levels/'+str(i) + '-' + str(i+1))
	p = p+1



x=1



