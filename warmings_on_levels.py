
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
from scipy.stats.stats import pearsonr





alt_constraint = iris.Constraint(Pressure = [30,50,70,100])


def callback(cube, field, filename): 
    # Remove the history attribute. 
    del cube.attributes['history']
    del cube.attributes['valid_max']
    del cube.attributes['valid_min']

fname = '../long_runs/CMIP_data/U/N96/daily/u-ar766_U_1850-2350.nc'
cubes = iris.load(fname, alt_constraint & lat_constraint, callback=callback)
U_N96 = cubes.concatenate_cube()
U_N96 = add_times(U_N96, 'time')



con1 = iris.Constraint(month = ['Nov', 'Dec'])
con2 = iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr'])
p = 0
for i in years_N96:
	U = []
	print i
	dummy2 = iris.cube.CubeList([U_N96_10.extract(iris.Constraint(year = i) & con1) ,U_N96_10.extract(iris.Constraint(year = i+1) & con2)])
	U.append(dummy2.concatenate_cube())
	plt.plot(np.arange(180), U[0].data)
	plt.plot([0,150], [0,0])
	plt.xticks([0,30,60,90,120,150], ['01/nov', '01/dec', '01/jan', '01/feb', '01/mar', '01/apr']
	plt.title(str(i) + '/' + str(i+1) + '   DJF mean U = ' + str(U_N96_10djf[p].data) + ' number of Warmings = ' + str(year_countN96[p]))
	plt.show()
	p = p+1
