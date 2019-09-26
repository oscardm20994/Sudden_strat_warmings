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
end_yearUKESM = 2559
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
fname = '../../long_runs/CMIP_data/U/N96/daily/u-ar766_U_zonmean10hPa_1850-2349.nc'
cubes = iris.load(fname, lat_constraint)#, callback=callback)
U_N96 = cubes[0]
U_N96 = add_times(U_N96, 'time')
U_N96 = U_N96[:,0,0]

fname = '../../long_runs/CMIP_data/U/N216/daily/u-aq281_U10hPazonmean*'
cubes = iris.load(fname, lat_constraint)#, callback = callback)
#cubes[-1] = cubes[-1].extract(iris.Constraint(Pressure = 10.))
U_N216 = cubes[0]
U_N216 = add_times(U_N216, 'time')
U_N216 = U_N216[:,0,0]

#UKESM U field
fname = '/gws/nopw/j04/aopp/strat_clim_data/CMIP6/UKESM1/picontrol/u_pl/daily/UKESM1*'
cubes = iris.load(fname, lat_constraint, callback=callback)
del cubes[3].attributes['nco_openmp_thread_number']
del cubes[3].attributes['NCO']
U_UKESM = cubes.concatenate_cube()
U_UKESM = add_times(U_UKESM, 't')
U_UKESM = U_UKESM[:,0,0]

#ERA reanlysis data
fname = '../../ERA_data/U/U_ERA_10hPa_daily.nc'
U_ERA_cubes = iris.load(fname, lat_constraint)
U_ERA = U_ERA_cubes[0]

def callback(cube, field, filename): 
    # Remove the history attribute. 
    del cube.attributes['history']
#ERA 40 data
fname = '../../ERA_data/U/U_ERA40/*.nc'
cubes = iris.load(fname, callback=callback)
del cubes[27]
U_ERA40 = cubes.concatenate_cube()
U_ERA40 = add_times(U_ERA40, 'time')

#ERA data from 2015
fname = '../../ERA_data/U/U_ERA_interim2/*.nc'
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


################################################################SSW counting and histograms begins here

months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']

#GC3
SSW_N216, SSW_year_N216, indexi216, indexj216 = SSW_counter(U_N216, start_yearN216, end_yearN216, 0)
SSW_ERA, SSW_year_ERA, indexiERA, indexjERA = SSW_counter(U_strat_ERA, start_yearERA, end_yearERA, 0)
SSW_ERA40, SSW_year_ERA40, indexiERA40, indexjERA40 = SSW_counter(U_strat_ERA40, start_yearERA40, end_yearERA40, 0)
SSW_ERA2, SSW_year_ERA2, indexiERA2, indexjERA2 = SSW_counter(U_strat_ERA2, start_yearERA2, end_yearERA2, 0)
SSW_UKESM, SSW_year_UKESM, indexiUK, indexjUK, time_list_UKESM = SSW_counter(U_UKESM, start_yearUKESM, end_yearUKESM, 0)
SSW_N96, SSW_year_N96, indexi96, indexj96 = SSW_counter_SP(U_N96, start_yearN96, end_yearN96, 0)

np.savetxt('event_list_UKESM.csv', np.array(time_list_UKESM), delimiter=",")




SSW_ERA.insert(0,'Feb')
SSW_ERA.insert(0,'Feb')

SSW_year_ERA.insert(0,'1980')
SSW_year_ERA.insert(0,'1979')

SSW_UKESM = ['Jan', 'Dec', 'Feb', 'Mar', 'Dec', 'Dec', 'Jan', 'Dec', 'Nov', 'Feb', 'Mar', 'Dec', 'Mar', 'Mar', 'Dec', 'Feb', 'Nov', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Jan', 'Nov', 'Feb', 'Nov', 'Mar', 'Mar', 'Dec', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Nov', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Nov', 'Feb', 'Feb', 'Mar', 'Mar', 'Mar', 'Mar', 'Feb', 'Jan', 'Nov', 'Feb', 'Feb', 'Nov', 'Mar', 'Feb', 'Jan', 'Jan', 'Dec', 'Jan', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Nov', 'Mar', 'Mar', 'Jan', 'Feb', 'Dec', 'Jan', 'Nov', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Feb', 'Mar', 'Dec', 'Mar', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Nov', 'Jan', 'Jan', 'Feb', 'Feb', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Feb', 'Mar', 'Jan', 'Feb', 'Nov', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Dec', 'Dec', 'Mar', 'Nov', 'Feb', 'Dec', 'Jan', 'Feb', 'Nov', 'Jan', 'Feb', 'Feb', 'Mar', 'Mar', 'Jan', 'Mar', 'Dec', 'Feb', 'Dec', 'Mar', 'Feb', 'Jan', 'Feb', 'Nov', 'Feb', 'Mar', 'Mar', 'Nov', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Mar', 'Jan', 'Feb', 'Mar', 'Jan', 'Dec', 'Feb', 'Mar', 'Feb', 'Jan', 'Jan', 'Nov', 'Jan', 'Mar', 'Jan', 'Dec', 'Mar', 'Feb', 'Mar', 'Feb', 'Nov', 'Jan', 'Jan', 'Mar', 'Feb', 'Jan', 'Nov', 'Feb', 'Jan', 'Nov', 'Nov', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Mar', 'Mar', 'Nov', 'Jan', 'Jan', 'Nov', 'Mar', 'Jan', 'Feb', 'Nov', 'Mar']

SSW_year_UKESM = [2061, 2063, 2063, 2063, 2064, 2067, 2067, 2069, 2072, 2076, 2079, 2082, 2083, 2085, 2089, 2092, 2093, 2094, 2095, 2096, 2100, 2101, 2107, 2110, 2111, 2111, 2113, 2114, 2115, 2116, 2118, 2121, 2121, 2122, 2122, 2123, 2126, 2129, 2130, 2131, 2131, 2132, 2135, 2135, 2140, 2142, 2143, 2144, 2145, 2146, 2149, 2150, 2151, 2152, 2153, 2154, 2156, 2159, 2162, 2164, 2165, 2169, 2169, 2170, 2171, 2173, 2174, 2175, 2179, 2179, 2179, 2180, 2182, 2182, 2183, 2187, 2191, 2192, 2194, 2197, 2201, 2205, 2209, 2209, 2212, 2213, 2214, 2215, 2217, 2218, 2219, 2220, 2222, 2223, 2227, 2228, 2228, 2231, 2234, 2238, 2245, 2246, 2246, 2252, 2254, 2257, 2258, 2259, 2259, 2260, 2260, 2265, 2266, 2267, 2269, 2270, 2270, 2272, 2277, 2277, 2280, 2280, 2280, 2282, 2283, 2285, 2285, 2286, 2288, 2289, 2291, 2293, 2294, 2294, 2297, 2297, 2300, 2312, 2313, 2315, 2315, 2316, 2321, 2322, 2325, 2329, 2329, 2330, 2332, 2335, 2335, 2337, 2339, 2341, 2341, 2343, 2346, 2346, 2347, 2349, 2352, 2354, 2356, 2356, 2360, 2362, 2365, 2368, 2369, 2374, 2376, 2376, 2379, 2379, 2380, 2381, 2392, 2393, 2393, 2395, 2395, 2398, 2404, 2404, 2405, 2415, 2416, 2418, 2421, 2422, 2425, 2425, 2426, 2428, 2429, 2430, 2439, 2441, 2443, 2444, 2446, 2447, 2448, 2450, 2450, 2456, 2458] 

SSW_N96 = ['Feb', 'Jan', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Mar', 'Dec', 'Mar', 'Feb', 'Mar', 'Jan', 'Jan', 'Nov', 'Nov', 'Feb', 'Feb', 'Feb', 'Feb', 'Feb', 'Dec', 'Dec', 'Nov', 'Feb', 'Mar', 'Jan', 'Mar', 'Nov', 'Mar', 'Nov', 'Nov', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Jan', 'Mar', 'Feb', 'Mar', 'Mar', 'Feb', 'Jan', 'Mar', 'Feb', 'Jan', 'Jan', 'Feb', 'Mar', 'Feb', 'Dec', 'Jan', 'Feb', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Nov', 'Feb', 'Jan', 'Dec', 'Jan', 'Jan', 'Feb', 'Mar', 'Jan', 'Feb', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Nov', 'Feb', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Feb', 'Mar', 'Feb', 'Dec', 'Nov', 'Feb', 'Jan', 'Jan', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Mar', 'Nov', 'Jan', 'Jan', 'Jan', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Nov', 'Dec', 'Feb', 'Jan', 'Mar', 'Nov', 'Dec', 'Feb', 'Feb', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Dec', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Feb', 'Mar', 'Dec', 'Jan', 'Mar', 'Mar', 'Nov', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Nov', 'Feb', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Mar', 'Nov', 'Feb', 'Jan', 'Nov', 'Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Mar', 'Feb', 'Mar', 'Jan', 'Feb', 'Feb', 'Nov', 'Dec', 'Feb', 'Mar', 'Jan', 'Dec', 'Feb', 'Dec', 'Mar', 'Feb', 'Mar', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Mar', 'Dec', 'Jan', 'Feb', 'Feb', 'Dec', 'Jan', 'Nov', 'Feb', 'Dec', 'Jan', 'Feb', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Dec', 'Feb', 'Mar', 'Dec', 'Mar', 'Jan', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Feb', 'Nov', 'Feb', 'Nov', 'Jan', 'Nov', 'Feb', 'Nov', 'Mar', 'Feb', 'Jan', 'Dec', 'Jan', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Nov', 'Mar', 'Mar', 'Dec', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Nov', 'Nov', 'Jan', 'Nov', 'Jan', 'Mar', 'Dec', 'Dec', 'Jan', 'Jan', 'Jan', 'Jan', 'Mar', 'Feb', 'Dec', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Dec', 'Dec', 'Feb', 'Mar', 'Mar', 'Jan', 'Jan', 'Nov', 'Dec', 'Nov', 'Jan', 'Feb', 'Mar', 'Dec', 'Jan', 'Jan', 'Mar', 'Jan', 'Nov', 'Nov', 'Feb', 'Mar', 'Mar', 'Feb', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Jan', 'Mar', 'Mar', 'Mar', 'Nov', 'Feb', 'Feb', 'Dec', 'Mar', 'Mar', 'Nov', 'Dec', 'Mar', 'Feb', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Jan', 'Jan', 'Feb', 'Dec', 'Mar', 'Jan', 'Feb', 'Nov', 'Nov', 'Mar', 'Jan', 'Mar', 'Dec', 'Feb']

SSW_year_N96 = [1852, 1854, 1854, 1855, 1855, 1856, 1857, 1857, 1862, 1862, 1863, 1864, 1865, 1865, 1868, 1869, 1870, 1872, 1873, 1873, 1874, 1875, 1876, 1878, 1879, 1880, 1880, 1881, 1882, 1883, 1885, 1886, 1886, 1887, 1888, 1888, 1889, 1891, 1892, 1894, 1895, 1896, 1896, 1897, 1901, 1902, 1903, 1904, 1904, 1905, 1906, 1907, 1908, 1909, 1911, 1912, 1913, 1915, 1918, 1918, 1920, 1920, 1922, 1922, 1923, 1924, 1926, 1927, 1929, 1929, 1930, 1933, 1935, 1938, 1939, 1941, 1941, 1943, 1944, 1946, 1946, 1947, 1948, 1948, 1949, 1950, 1950, 1951, 1952, 1954, 1954, 1955, 1955, 1957, 1958, 1958, 1960, 1960, 1961, 1964, 1966, 1967, 1968, 1970, 1973, 1974, 1976, 1976, 1977, 1979, 1980, 1982, 1982, 1987, 1988, 1989, 1990, 1991, 1991, 1991, 1992, 1992, 1994, 1994, 1996, 1997, 2000, 2001, 2001, 2001, 2002, 2004, 2005, 2008, 2009, 2010, 2012, 2013, 2013, 2014, 2015, 2017, 2018, 2019, 2019, 2024, 2024, 2025, 2028, 2029, 2030, 2030, 2034, 2039, 2040, 2043, 2044, 2045, 2046, 2046, 2049, 2052, 2054, 2057, 2057, 2058, 2059, 2060, 2061, 2066, 2067, 2070, 2071, 2072, 2073, 2073, 2074, 2075, 2075, 2076, 2077, 2077, 2081, 2085, 2087, 2087, 2088, 2090, 2091, 2093, 2094, 2094, 2097, 2097, 2099, 2100, 2100, 2103, 2103, 2105, 2105, 2106, 2107, 2107, 2107, 2108, 2108, 2116, 2116, 2118, 2119, 2122, 2122, 2123, 2123, 2123, 2124, 2125, 2125, 2125, 2126, 2126, 2128, 2130, 2131, 2131, 2132, 2132, 2134, 2137, 2138, 2145, 2146, 2146, 2146, 2149, 2150, 2150, 2151, 2153, 2155, 2156, 2158, 2160, 2160, 2161, 2161, 2163, 2164, 2169, 2169, 2171, 2173, 2176, 2176, 2178, 2179, 2179, 2180, 2181, 2181, 2181, 2182, 2183, 2186, 2187, 2188, 2188, 2189, 2191, 2192, 2193, 2194, 2194, 2198, 2200, 2201, 2201, 2203, 2204, 2206, 2209, 2211, 2215, 2215, 2215, 2216, 2217, 2218, 2220, 2220, 2220, 2221, 2222, 2223, 2224, 2224, 2225, 2226, 2228, 2229, 2230, 2230, 2231, 2233, 2238, 2240, 2244, 2245, 2246, 2246, 2248, 2249, 2250, 2250, 2252, 2252, 2253, 2255, 2255, 2256, 2257, 2258, 2259, 2259, 2263, 2264, 2264, 2266, 2266, 2268, 2268, 2269, 2272, 2274, 2275, 2275, 2275, 2276, 2277, 2278, 2279, 2283, 2284, 2287, 2291, 2292, 2294, 2295, 2295, 2296, 2302, 2304, 2306, 2307, 2307, 2308, 2308, 2309, 2310, 2310, 2313, 2313, 2314, 2316, 2319, 2320, 2320, 2322, 2322, 2326, 2327, 2330, 2331, 2332, 2337, 2337, 2339, 2340, 2341, 2342, 2343, 2344, 2345, 2345, 2346, 2349]

SSW_N216 = ['Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Dec', 'Feb', 'Jan', 'Feb', 'Feb', 'Nov', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan', 'Mar', 'Mar', 'Mar', 'Feb', 'Nov', 'Nov', 'Jan', 'Jan', 'Mar', 'Feb', 'Feb', 'Mar', 'Feb', 'Jan', 'Mar', 'Jan', 'Mar', 'Feb', 'Dec', 'Feb', 'Mar', 'Jan', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Mar', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Nov', 'Mar', 'Mar', 'Jan', 'Mar', 'Feb', 'Feb', 'Jan', 'Jan', 'Feb', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Nov', 'Feb', 'Jan', 'Jan', 'Feb', 'Dec', 'Feb', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Mar', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Jan', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Dec', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Nov', 'Jan', 'Feb', 'Feb', 'Mar', 'Dec', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Feb', 'Jan', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb']

SSW_year_N216 = [1852, 1855, 1857, 1858, 1860, 1860, 1861, 1862, 1862, 1863, 1863, 1864, 1865, 1866, 1870, 1872, 1875, 1880, 1884, 1885, 1889, 1891, 1891, 1892, 1893, 1893, 1897, 1898, 1898, 1899, 1899, 1900, 1902, 1909, 1911, 1915, 1916, 1918, 1919, 1922, 1922, 1923, 1924, 1926, 1927, 1928, 1929, 1929, 1930, 1934, 1934, 1935, 1936, 1937, 1937, 1939, 1939, 1942, 1942, 1943, 1944, 1944, 1946, 1947, 1948, 1951, 1952, 1953, 1954, 1957, 1963, 1965, 1966, 1967, 1967, 1967, 1969, 1970, 1970, 1971, 1973, 1974, 1977, 1980, 1980, 1981, 1982, 1985, 1988, 1989, 1990, 1992, 1993, 1995, 1998, 2002, 2003, 2010, 2012, 2013, 2015, 2017, 2017, 2020, 2020, 2020, 2021, 2022, 2023, 2023, 2023, 2025, 2029, 2029, 2030, 2032, 2035, 2038, 2039, 2040, 2041, 2042, 2048, 2048, 2048, 2049, 2050, 2051, 2053, 2053, 2054, 2056, 2058, 2059, 2059, 2060, 2061, 2064, 2070, 2071, 2072, 2073, 2075, 2077, 2078, 2079, 2081, 2081, 2082, 2083, 2084, 2084, 2085, 2087, 2087, 2089, 2090, 2091, 2093, 2094, 2096, 2099, 2099, 2101, 2102, 2103, 2104, 2105, 2107, 2108, 2109, 2111, 2113, 2114, 2115, 2118, 2121, 2122, 2123, 2123, 2124, 2127, 2127, 2129, 2132, 2137, 2138, 2139, 2140, 2141, 2143, 2145, 2145, 2149] 

SSW_ERA = ['Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Jan', 'Jan', 'Dec', 'Mar', 'Feb', 'Dec', 'Feb', 'Mar', 'Feb', 'Dec', 'Jan', 'Jan', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan']

SSW_year_ERA = [1979, 1980, 1981, 1981, 1984, 1985, 1987, 1987, 1988, 1989, 1998, 1999, 2000, 2001, 2001, 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2010, 2013]

SSW_ERA40 = ['Jan', 'Dec', 'Feb', 'Jan', 'Nov', 'Mar', 'Jan', 'Jan', 'Mar', 'Jan', 'Jan']
SSW_year_ERA40 = [1963, 1965, 1966, 1968, 1968, 1969, 1970, 1971, 1971, 1973, 1977]

SSW_ERA2 = ['Feb', 'Feb', 'Mar']
SSW_year_ERA2 = [2017, 2018, 2018]

SSW_ERA = SSW_ERA40 + SSW_ERA + SSW_ERA2
SSW_year_ERA = SSW_year_ERA40 + SSW_year_ERA + SSW_year_ERA2

#count SSWs for each month
SSW_N96_counter = Counter(SSW_N96)
SSW_N216_counter = Counter(SSW_N216)
SSW_ERA_counter = Counter(SSW_ERA)
SSW_UKESM_counter = Counter(SSW_UKESM)


########## Warmings per decade by taking average and STDEV of counts per decade
per_decadeERA, ErERA, running_ERA, year_countERA = DJF_warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA)
per_decadeN96, ErN96, running_N96, year_countN96 = DJF_warmings_per_decade(SSW_N96, SSW_year_N96, years_N96)
per_decadeN216, ErN216, running_N216, year_countN216 = DJF_warmings_per_decade(SSW_N216, SSW_year_N216, years_N216)
per_decadeUKESM, ErUKESM, running_UKESM, year_countUKESM = DJF_warmings_per_decade(SSW_UKESM, SSW_year_UKESM, years_UKESM)

per_decade = [per_decadeERA, per_decadeN96, per_decadeN216]
Er = [ErERA, ErN96, ErN216]


########## Warmings per decade for each month
start_yearERAc = 1959 
end_yearERAc = 2018
years_ERAc = np.arange(start_yearERAc, end_yearERAc, 1)


N96, N216, ERA, UKESM, ErN96, ErN216, ErERA, ErUKESM = [], [], [], [], [], [], [], []
for i in range(len(months)):
	N96.append(mon_warmings_per_decade(SSW_N96, SSW_year_N96, years_N96, months[i])[0])
	UKESM.append(mon_warmings_per_decade(SSW_UKESM, SSW_year_UKESM, years_UKESM, months[i])[0])
	N216.append(mon_warmings_per_decade(SSW_N216, SSW_year_N216, years_N216, months[i])[0])
	ERA.append(mon_warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA, months[i])[0])
	ErN96.append(mon_warmings_per_decade(SSW_N96, SSW_year_N96, years_N96, months[i])[1])
	ErN216.append(mon_warmings_per_decade(SSW_N216, SSW_year_N216, years_N216, months[i])[1])
	ErERA.append(mon_warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA, months[i])[1])
	ErUKESM.append(mon_warmings_per_decade(SSW_UKESM, SSW_year_UKESM, years_UKESM, months[i])[1])


import pickle


picklefile='Histogramvariables.pickle'
# Output variables used in the plot
pfile = open(picklefile,'wb')
pickle.dump([N96,N216,ERA,ErN96,ErN216,ErERA],pfile)
pfile.close()


index = np.arange(5)
fig, ax = plt.subplots(figsize=(18,12))
bar_width = 0.15
warm1 = ax.bar(index, N96, bar_width, yerr=ErN96,  capsize = 10, ecolor = 'black', color = 'r', linewidth = 2.0)
warm2 = ax.bar(index + bar_width, N216, bar_width, yerr=ErN216, capsize = 10, ecolor = 'black', color = 'g', linewidth = 2.0)
warm4 = ax.bar(index + 3*bar_width, UKESM, bar_width, yerr=ErUKESM, capsize = 10, ecolor = 'black', color = 'b', linewidth = 2.0)
warm3 = ax.bar(index + 2*bar_width, ERA, bar_width, yerr=ErERA, capsize = 10, ecolor = 'black', color = 'cyan', linewidth = 2.0)
ax.set_xticks(index + 2*bar_width)
ax.set_xticklabels((months), fontsize = 20)
ax.legend((warm1[0],warm2[0],warm3[0], warm4[0]),('HadgemGC3.1 N96','HadgemGC3.1 N216','ERA','UKESM1'), fontsize = 20, loc = 2)
plt.title('Mean SSWs per season', fontsize = 30)
plt.ylabel('SSWs per season', fontsize = 20)
plt.yticks(np.arange(6), fontsize = 20)
plt.xlim(-0.25,5)
plt.ylim(0,0.4) 
plt.show()
fig.savefig('../figures/SSw_perdecade_months.png')


month = [' ','Jul',' ', 'Aug',' ', 'Sep',' ', 'Oct',' ', 'Nov' ,' ', 'Dec',' ', 'Jan',' ', 'Feb',' ', 'Mar',' ', 'Apr',' ', 'May',' ', 'Jun']

#take day means
N96_daymeans1 = U_N96.aggregated_by('day_number', iris.analysis.MEAN)
N216_daymeans1 = U_N216.aggregated_by('day_number', iris.analysis.MEAN)
UKESM_daymeans1 = U_UKESM.aggregated_by('day_number', iris.analysis.MEAN)

ERA_daymeans1 = U_strat_ERA.aggregated_by('day_number', iris.analysis.MEAN)
U_strat_ERA402 = iris.cube.CubeList([U_strat_ERA40, U_strat_ERA2]).concatenate_cube()
ERA402_daymeans1 = U_strat_ERA402.aggregated_by('day_number', iris.analysis.MEAN)

#day standard deviations
N96_daystd1 = U_N96.aggregated_by('day_number', iris.analysis.STD_DEV)
N216_daystd1 = U_N216.aggregated_by('day_number', iris.analysis.STD_DEV)
UKESM_daystd1 = U_UKESM.aggregated_by('day_number', iris.analysis.STD_DEV)

ERA_daystd1 = U_strat_ERA.aggregated_by('day_number', iris.analysis.STD_DEV)
ERA402_daystd1 = U_strat_ERA402.aggregated_by('day_number', iris.analysis.STD_DEV)

N96_daymean = np.append(N96_daymeans1[180:360].data, N96_daymeans1[0:180].data)
N216_daymean = np.append(N216_daymeans1[180:360].data, N216_daymeans1[0:180].data)
UKESM_daymean = np.append(UKESM_daymeans1[180:360].data, UKESM_daymeans1[0:180].data)
ERA_daymean = np.append(ERA_daymeans1[180:360].data, ERA_daymeans1[0:180].data)
ERA402_daymean = np.append(ERA402_daymeans1[180:360].data, ERA402_daymeans1[0:180].data)

ERA_daymeanc = np.mean([ERA_daymean, ERA402_daymean], axis = 0)

N96_daystd = np.append(N96_daystd1[180:360].data, N96_daystd1[0:180].data)
N216_daystd = np.append(N216_daystd1[180:360].data, N216_daystd1[0:180].data)
UKESM_daystd = np.append(UKESM_daystd1[180:360].data, UKESM_daystd1[0:180].data)
ERA_daystd = np.append(ERA_daystd1[180:360].data, ERA_daystd1[0:180].data)
ERA402_daystd = np.append(ERA402_daystd1[180:360].data, ERA402_daystd1[0:180].data)

ERA_daystdc = np.mean([ERA_daystd, ERA402_daystd], axis = 0)
ERA_daystdc = 0.5*np.sqrt(ERA_daystd**2 + ERA402_daystd**2)
time1 = np.linspace(0,12,len(N96_daymean))
time2 = np.linspace(0,12,len(ERA_daymean))

import pickle
picklefile='lineplotvariables.pickle'
pfile = open(picklefile,'wb')
pickle.dump([N96_daymean,N216_daymean,ERA_daymeanc,N96_daystd,N216_daystd,ERA_daystdc],pfile)
pfile.close()


fig, ax = plt.subplots(dpi=300)
plt.plot([0,12], [0,0], linewidth = 0.75)
plt.plot(time1, N96_daymean, color= 'r', label = 'GC3 N96', linewidth = 0.75)
plt.plot(time1, UKESM_daymean, color= 'b', label = 'UKESM ', linewidth = 0.75)
#plt.plot(time1, N216_daymean, color= 'g', label = 'GC3 N216 ', linewidth = 0.75)
plt.plot(time2, ERA_daymean, color= 'black', label = 'ERA', linewidth = 0.9)
#plt.plot(time1, N216_daymean + N216_daystd, linestyle = '--', color = 'g', linewidth = 0.75)
#plt.plot(time1, N216_daymean - N216_daystd, linestyle = '--', color = 'g', linewidth = 0.75)
plt.plot(time1, N96_daymean + N96_daystd, linestyle = '--', color = 'r', linewidth = 0.75)
plt.plot(time1, N96_daymean - N96_daystd, linestyle = '--', color = 'r', linewidth = 0.75)
plt.plot(time1, UKESM_daymean + UKESM_daystd, linestyle = '--', color = 'b', linewidth = 0.75)
plt.plot(time1, UKESM_daymean - UKESM_daystd, linestyle = '--', color = 'b', linewidth = 0.75)
plt.plot(time2, ERA_daymean + ERA_daystd, linestyle = '--', color = 'black', linewidth = 0.75)
plt.plot(time2, ERA_daymean - ERA_daystd, linestyle = '--', color = 'black', linewidth = 0.75)

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
plt.ylabel('$ZMZW (ms^{-1})$')
plt.legend()
plt.xlim(0,12)
plt.title('60N 10hPa ZMZW climatology')
plt.show()
fig.savefig('../../figures/60s_10hPa_resoloutions_for_ACSIS.png')


# t test
scipy.stats.ttest_ind(a, b, axis=0, equal_var=False)



#running 10 year count
fig = plt.figure(figsize=(18, 9))
plt.plot(years_N96[0:-10], running_N96, label = 'N96')
#plt.plot(years_N216[0:-10], running_N216, label = 'N216')
plt.title('DJF SSWs per decade, HadgemGC3.1 N96 resoloution', fontsize = 20)
plt.xlabel('Decade start', fontsize = 20)
plt.ylabel('SSWs per decade', fontsize = 20)
plt.xticks(fontsize = 13)
plt.yticks(fontsize = 13)
plt.xlim(years_N96[0],years_N96[-10])
plt.show()
fig.savefig('../figures/running_SSW_count_N96.png')





year_count = np.empty(0)
year_count_nov = np.empty(0)
year_random = np.zeros((len(years_GC3), 5))
SSWs_GC3 = [SSW_GC3_counter['Nov'],SSW_GC3_counter['Dec'],SSW_GC3_counter['Jan'],SSW_GC3_counter['Feb'],SSW_GC3_counter['Mar']]

#GC3 SSW count seperated by year
fig = plt.figure(figsize=(18, 9))
plt.plot(years_UKESM[0:-10], running_UKESM)
#plt.xlim(years_UKESM[0], years_UKESM[-10])
plt.yticks([0,1,2,3,4,5])
plt.title('DJF SSWs per year')
plt.ylabel('number of SSWs')
plt.xlabel('year')
plt.show()
fig.savefig('../figures/year_count_N96.png')

fig = plt.figure(figsize=(18, 9))
plt.bar(years_UKESM, year_countUKESM)
#plt.bar(years_N216, year_countN96)
plt.ylim(-0.1,5)
plt.yticks([0,1,2,3,4,5])
plt.title('DJF SSWs per year')
plt.ylabel('number of SSWs')
plt.xlabel('year')
plt.show()
fig.savefig('../figures/DJF_SSWs_per_year_N216.png')


















