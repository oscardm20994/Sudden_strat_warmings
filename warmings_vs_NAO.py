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
from eofs.iris import Eof
from eofs.standard import Eof as Eof2

################################################READING IN DATA AND INITIALISING ARRAYS

#start and end years for data
start_yearN96 = 1850
end_yearN96 = 2351
years_N96 = np.arange(start_yearN96, end_yearN96, 1)

### EOF region
def correct_lat_EOF(cell):
	return 19.9 < cell < 80.1

lat_constraint_EOF = iris.Constraint(latitude = correct_lat_EOF)

def correct_lon_EOF(cell):
    return    269.9 < cell < 360 or 0 < cell < 40.1

lon_constraint_EOF = iris.Constraint(longitude = correct_lon_EOF)

### jet lat and lon 
def correct_lat_jet(cell):
	return 14.9 < cell < 75.1

lat_constraint_jet = iris.Constraint(latitude = correct_lat_jet)

def correct_lon_jet(cell):
    return    299.9 < cell < 360.1

#station region
lon_constraint_jet = iris.Constraint(longitude = correct_lon_jet)
lon_lisbon = iris.Constraint(longitude = [351.5625 , 350.])
lat_lisbon = iris.Constraint(latitude = [38.125 , 38.])
lat_ice = iris.Constraint(longitude = [338.4375 , 338.])
lon_ice = iris.Constraint(latitude = [64.375 , 64.])

##################################### define functions
def callback(cube, field, filename): 
    # Remove the history attribute. 
    del cube.attributes['history']
    del cube.attributes['valid_max']
    del cube.attributes['valid_min']

def NAO_stations(SLP, time):
	SLP_lisb = SLP.extract(lat_lisbon & lon_lisbon)
	SLP_ice = SLP.extract(lat_ice & lon_ice)
	SLP_icedjf = SLP_ice.extract(iris.Constraint(clim_season = 'djf'))
	SLP_icedjf = SLP_icedjf.aggregated_by('season_year', iris.analysis.MEAN)
	SLP_lisbdjf = SLP_lisb.extract(iris.Constraint(clim_season = 'djf'))
	SLP_lisbdjf = SLP_lisbdjf.aggregated_by('season_year', iris.analysis.MEAN)
	print SLP_icedjf
	meanice = SLP_icedjf.collapsed(time, iris.analysis.MEAN)
	stdevice = SLP_icedjf.collapsed(time, iris.analysis.STD_DEV)
	meanlisb = SLP_lisbdjf.collapsed(time, iris.analysis.MEAN)
	stdevlisb = SLP_lisbdjf.collapsed(time, iris.analysis.STD_DEV)
	SLP_icedjf = (SLP_icedjf - meanice)/stdevice
	SLP_lisbdjf = (SLP_lisbdjf - meanlisb)/stdevice
	NAO_stations = SLP_lisbdjf - SLP_icedjf
	return NAO_stations

def NAO_EOFs(SLP, neofs, time, months):
	SLP_NAO = SLP.extract(lat_constraint_EOF & lon_constraint_EOF)
	SLP_djf = SLP_NAO.extract(iris.Constraint(month = months))
	SLP_djf = SLP_djf.aggregated_by('year', iris.analysis.MEAN)
	mean = SLP_djf.collapsed(time, iris.analysis.MEAN)
	stdev = SLP_djf.collapsed(time, iris.analysis.STD_DEV)
	SLP_djf = (SLP_djf - mean)/stdev
	solver = Eof(SLP_djf, weights = 'coslat')
	eof = solver.eofs(neofs=neofs)
	NAO_djf = solver.pcs(npcs=neofs, pcscaling=1)[:,0]
	return NAO_djf, eof

def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.
    Args:
    window: int
        The length of the filter window.
    cutoff: float
        The cutoff frequency in inverse time steps.
    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]

#################load data
#GC3 U fields
fname = '../long_runs/CMIP_data/U/N96/daily/u-ar766_U_z*.nc'
cubes = iris.load(fname, lat_constraint)
U_N96_10 = cubes[0][:,0,0]
U_N96_10 = add_times(U_N96_10, 'time')
SSW_N96, SSW_year_N96, indexi96, indexj96 = SSW_counter(U_N96_10, years_N96[0], years_N96[-1], 0)

#U fields for jet strength
fname = '../long_runs/CMIP_data/U/N96/daily/u-ar766_U_1850-2350.nc'
cubes = iris.load(fname, lat_constraint_jet & lon_constraint_jet)
U_N96_jet = cubes[0][:,0:3,:,:]
U_N96_jet = U_N96_jet.collapsed(['longitude', 'air_pressure'], iris.analysis.MEAN)

#GC3 SLP field
fname = '../long_runs/CMIP_data/SLP/N96/monthly/*.nc'
SLP_N96 = iris.load(fname)[0]
SLP_N96 = add_times(SLP_N96, 't')
SLP_N96 = SLP_N96[:,0,:,:]
SLP_GC3 = SLP_N96
print SLP_N96

#ERA SLP data
fname = '../ERA_data/SLP_monthly/*'
cubes = iris.load(fname, callback=callback)
SLP_ERA = cubes.concatenate_cube()
SLP_ERA = add_times(SLP_ERA, 'time')




###############################################ANALYSIS STARTS HERE ####################
##NAO calculation using EOF and stations
NAO = NAO_EOFs(SLP_N96, 1, 't')[0]
NAO_sta = NAO_stations(SLP_N96, 't')

NAO_ERA = NAO_EOFs(SLP_ERA, 1, 'time')[0]
NAO_ERA_sta = NAO_stations(SLP_ERA, 'time')

NAO_obs = np.loadtxt('NAO_obs.txt', dtype = float)
NAO_obsdjf = np.array([NAO_obs[:,1],NAO_obs[:,2],NAO_obs[:,-2]])
NAO_obsmn = np.mean(NAO_obsdjf, axis = 0)

fig = plt.figure(figsize=(18, 9)) 
plt.bar(NAO_ERA.coord('season_year').points, NAO_ERA.data, label = 'PC derived')
plt.plot(NAO_ERA.coord('season_year').points, NAO_obsmn, label = 'observations', color = 'g')
plt.plot(NAO_ERA.coord('season_year').points, NAO_ERA_sta.data, color = 'r', label = 'station derived')
plt.title('NAO index from ERA interim dataset')
plt.legend()
plt.show()






fig = plt.figure(figsize=(18, 9)) 
qplt.contourf(eofN96[0,:,:], cmap = 'viridis', levels = np.linspace(np.min(eofN96.data),np.max(eofN96.data), 20))
plt.gca().coastlines()
plt.title('1st EOF of sea level pressure, Hadgem GC3.1 pi control run, 1850-2349')
qplt.show()
fig.savefig('../figures/EOF1_SLP_N96.png')

#find jet strength
w = low_pass_weights(61, 0.1)
U_jet_filt = U_N96_jet.rolling_window('time', iris.analysis.MEAN, len(w), weights=w)
jet_speed = U_jet_filt.collapsed('latitude', iris.analysis.MAX)
jet_speed = add_times(jet_speed, 'time')
jet_speed_nofilt = U_N96_jet.collapsed('latitude', iris.analysis.MEAN)
jet_speed_nofilt = add_times(jet_speed_nofilt, 'time')



########################


## djf aggregates for U
U_N96_10djf = U_N96_10.extract(iris.Constraint(clim_season = 'djf'))
U_N96_10djf = U_N96_10djf.aggregated_by('season_year', iris.analysis.MEAN)

jet_djf = jet_speed.extract(iris.Constraint(clim_season = 'djf'))
jet_djf = jet_djf.aggregated_by('season_year', iris.analysis.MEAN)

jet_djf_nofilt = jet_speed_nofilt.extract(iris.Constraint(clim_season = 'djf'))
jet_djf_nofilt = jet_djf_nofilt.aggregated_by('season_year', iris.analysis.MEAN)

fig = plt.figure(figsize=(9, 9)) 
plt.scatter(NAO_lagged, U_SSW)
plt.title('ZMZW 60N 10hPa in month of SSW vs 3 month mean NAO after SSW, r = ' + str(np.corrcoef(NAO_lagged, U_SSW)[0,1]))
plt.xlabel('3 month lag NAO')
plt.ylabel('U 10hPa in SSW month')
plt.show()
fig.savefig('../figures/NAO_vs_10hPaU_3month.png')

SSW_N96 = ['Jan', 'Feb', 'Jan', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Mar', 'Dec', 'Mar', 'Feb', 'Mar', 'Jan', 'Jan', 'Nov', 'Nov', 'Feb', 'Feb', 'Feb', 'Feb', 'Feb', 'Dec', 'Dec', 'Nov', 'Feb', 'Mar', 'Jan', 'Mar', 'Nov', 'Mar', 'Nov', 'Nov', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Jan', 'Mar', 'Feb', 'Mar', 'Mar', 'Feb', 'Jan', 'Mar', 'Feb', 'Jan', 'Jan', 'Feb', 'Mar', 'Feb', 'Dec', 'Jan', 'Feb', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Nov', 'Feb', 'Jan', 'Dec', 'Jan', 'Jan', 'Feb', 'Mar', 'Jan', 'Feb', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Nov', 'Feb', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Feb', 'Mar', 'Feb', 'Dec', 'Nov', 'Feb', 'Jan', 'Jan', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Mar', 'Nov', 'Jan', 'Jan', 'Jan', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Nov', 'Dec', 'Feb', 'Jan', 'Mar', 'Nov', 'Dec', 'Feb', 'Feb', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Dec', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Feb', 'Mar', 'Dec', 'Jan', 'Mar', 'Mar', 'Nov', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Nov', 'Feb', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Mar', 'Nov', 'Feb', 'Jan', 'Nov', 'Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Mar', 'Feb', 'Mar', 'Jan', 'Feb', 'Feb', 'Nov', 'Dec', 'Feb', 'Mar', 'Jan', 'Dec', 'Feb', 'Dec', 'Mar', 'Feb', 'Mar', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Mar', 'Dec', 'Jan', 'Feb', 'Feb', 'Dec', 'Jan', 'Nov', 'Feb', 'Dec', 'Jan', 'Feb', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Dec', 'Feb', 'Mar', 'Dec', 'Mar', 'Jan', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Feb', 'Nov', 'Feb', 'Nov', 'Jan', 'Nov', 'Feb', 'Nov', 'Mar', 'Feb', 'Jan', 'Dec', 'Jan', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Nov', 'Mar', 'Mar', 'Dec', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Nov', 'Nov', 'Jan', 'Nov', 'Jan', 'Mar', 'Dec', 'Dec', 'Jan', 'Jan', 'Jan', 'Jan', 'Mar', 'Feb', 'Dec', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Dec', 'Dec', 'Feb', 'Mar', 'Mar', 'Jan', 'Jan', 'Nov', 'Dec', 'Nov', 'Jan', 'Feb', 'Mar', 'Dec', 'Jan', 'Jan', 'Mar', 'Jan', 'Nov', 'Nov', 'Feb', 'Mar', 'Mar', 'Feb', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Jan', 'Mar', 'Mar', 'Mar', 'Nov', 'Feb', 'Feb', 'Dec', 'Mar', 'Mar', 'Nov', 'Dec', 'Mar', 'Feb', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Jan', 'Jan', 'Feb', 'Dec', 'Mar', 'Jan', 'Feb', 'Nov', 'Nov', 'Mar', 'Jan', 'Mar', 'Dec', 'Feb']
SSW_year_N96 = [1851, 1852, 1854, 1854, 1855, 1855, 1856, 1857, 1857, 1862, 1862, 1863, 1864, 1865, 1866, 1868, 1869, 1870, 1872, 1873, 1873, 1874, 1875, 1876, 1878, 1879, 1880, 1881, 1882, 1882, 1883, 1885, 1886, 1886, 1887, 1888, 1888, 1889, 1891, 1892, 1894, 1895, 1896, 1896, 1897, 1901, 1902, 1903, 1904, 1904, 1905, 1906, 1907, 1908, 1909, 1911, 1912, 1913, 1915, 1918, 1918, 1920, 1921, 1922, 1922, 1923, 1924, 1926, 1927, 1929, 1929, 1930, 1933, 1936, 1938, 1939, 1941, 1941, 1943, 1944, 1946, 1946, 1947, 1948, 1948, 1949, 1950, 1950, 1951, 1952, 1954, 1954, 1955, 1956, 1957, 1958, 1958, 1960, 1961, 1961, 1964, 1966, 1967, 1968, 1970, 1973, 1974, 1976, 1976, 1977, 1979, 1980, 1982, 1982, 1987, 1988, 1989, 1990, 1991, 1991, 1991, 1992, 1992, 1994, 1994, 1997, 1997, 2000, 2001, 2001, 2002, 2002, 2004, 2005, 2008, 2009, 2010, 2012, 2013, 2013, 2015, 2015, 2017, 2018, 2019, 2019, 2024, 2024, 2026, 2028, 2029, 2030, 2030, 2034, 2039, 2040, 2043, 2044, 2045, 2046, 2046, 2049, 2052, 2054, 2057, 2057, 2058, 2059, 2060, 2061, 2066, 2067, 2070, 2071, 2072, 2073, 2073, 2074, 2075, 2075, 2076, 2077, 2077, 2081, 2085, 2087, 2088, 2088, 2090, 2091, 2093, 2094, 2094, 2097, 2097, 2100, 2100, 2100, 2103, 2104, 2105, 2106, 2106, 2107, 2107, 2107, 2108, 2108, 2116, 2116, 2118, 2119, 2122, 2122, 2123, 2123, 2124, 2124, 2125, 2125, 2126, 2126, 2126, 2128, 2131, 2131, 2131, 2132, 2132, 2134, 2137, 2138, 2146, 2146, 2146, 2147, 2149, 2150, 2150, 2151, 2153, 2155, 2156, 2158, 2160, 2160, 2161, 2161, 2163, 2164, 2169, 2169, 2171, 2173, 2176, 2177, 2178, 2179, 2179, 2180, 2181, 2181, 2181, 2182, 2183, 2186, 2187, 2188, 2188, 2189, 2191, 2193, 2193, 2194, 2194, 2198, 2200, 2201, 2201, 2203, 2204, 2206, 2209, 2211, 2215, 2215, 2216, 2217, 2217, 2218, 2220, 2220, 2220, 2221, 2223, 2223, 2224, 2224, 2225, 2226, 2228, 2229, 2230, 2230, 2232, 2234, 2238, 2240, 2244, 2245, 2246, 2246, 2249, 2249, 2250, 2250, 2252, 2253, 2253, 2255, 2255, 2256, 2257, 2258, 2259, 2259, 2263, 2264, 2264, 2266, 2266, 2268, 2268, 2269, 2272, 2274, 2275, 2275, 2275, 2276, 2277, 2278, 2279, 2283, 2284, 2287, 2292, 2292, 2294, 2295, 2296, 2296, 2302, 2304, 2306, 2307, 2307, 2308, 2308, 2309, 2310, 2310, 2313, 2314, 2314, 2316, 2319, 2320, 2320, 2322, 2322, 2326, 2327, 2330, 2331, 2332, 2337, 2338, 2339, 2340, 2341, 2342, 2343, 2344, 2345, 2345, 2347, 2349]
months = ['Nov' , 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct',]


months = 
NAO_lagged = np.empty(0)
U_SSW = np.empty(0)
for i in range(len(SSW_N96)):
	print i, SSW_N96[i]
	U = U_N96_10.extract(iris.Constraint(month = SSW_N96[i], year = SSW_year_N96[i]))
	U = U.collapsed('time', iris.analysis.MEAN)
	Us = np.append(Us, U)
	ind = months.index(SSW_N96[i])
	print ind
	NAO_lag = NAO_EOFs(SLP_N96, 1, 't', months[ind:ind+1])[0]
	NAO = NAO_lag.extract(iris.Constraint(year = SSW_year_N96[i]))
	NAOs = np.append(NAOs, NAO)

for i in range(388):
	NAO_lagged = np.append(NAO_lagged, NAOs[i])
	U_SSW = np.append(U_SSW, Us[i])




