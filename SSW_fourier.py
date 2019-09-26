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

#start and end years for data
start_yearN96 = 1850
end_yearN96 = 2349
years_N96 = np.arange(start_yearN96, end_yearN96, 1)

start_yearN216 = 1850
end_yearN216 = 2129
years_N216 = np.arange(start_yearN216, end_yearN216, 1)


SSW_N96 = ['Feb', 'Jan', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Mar', 'Dec', 'Mar', 'Feb', 'Mar', 'Jan', 'Jan', 'Nov', 'Nov', 'Feb', 'Feb', 'Feb', 'Feb', 'Feb', 'Dec', 'Dec', 'Nov', 'Feb', 'Mar', 'Jan', 'Mar', 'Nov', 'Mar', 'Nov', 'Nov', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Jan', 'Mar', 'Feb', 'Mar', 'Mar', 'Feb', 'Jan', 'Mar', 'Feb', 'Jan', 'Jan', 'Feb', 'Mar', 'Feb', 'Dec', 'Jan', 'Feb', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Nov', 'Feb', 'Jan', 'Dec', 'Jan', 'Jan', 'Feb', 'Mar', 'Jan', 'Feb', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Nov', 'Feb', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Feb', 'Mar', 'Feb', 'Dec', 'Nov', 'Feb', 'Jan', 'Jan', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Mar', 'Nov', 'Jan', 'Jan', 'Jan', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Nov', 'Dec', 'Feb', 'Jan', 'Mar', 'Nov', 'Dec', 'Feb', 'Feb', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Dec', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Feb', 'Mar', 'Dec', 'Jan', 'Mar', 'Mar', 'Nov', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Nov', 'Feb', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Mar', 'Nov', 'Feb', 'Jan', 'Nov', 'Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Mar', 'Feb', 'Mar', 'Jan', 'Feb', 'Feb', 'Nov', 'Dec', 'Feb', 'Mar', 'Jan', 'Dec', 'Feb', 'Dec', 'Mar', 'Feb', 'Mar', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Mar', 'Dec', 'Jan', 'Feb', 'Feb', 'Dec', 'Jan', 'Nov', 'Feb', 'Dec', 'Jan', 'Feb', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Dec', 'Feb', 'Mar', 'Dec', 'Mar', 'Jan', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Feb', 'Nov', 'Feb', 'Nov', 'Jan', 'Nov', 'Feb', 'Nov', 'Mar', 'Feb', 'Jan', 'Dec', 'Jan', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Nov', 'Mar', 'Mar', 'Dec', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Nov', 'Nov', 'Jan', 'Nov', 'Jan', 'Mar', 'Dec', 'Dec', 'Jan', 'Jan', 'Jan', 'Jan', 'Mar', 'Feb', 'Dec', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Dec', 'Dec', 'Feb', 'Mar', 'Mar', 'Jan', 'Jan', 'Nov', 'Dec', 'Nov', 'Jan', 'Feb', 'Mar', 'Dec', 'Jan', 'Jan', 'Mar', 'Jan', 'Nov', 'Nov', 'Feb', 'Mar', 'Mar', 'Feb', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Jan', 'Mar', 'Mar', 'Mar', 'Nov', 'Feb', 'Feb', 'Dec', 'Mar', 'Mar', 'Nov', 'Dec', 'Mar', 'Feb', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Jan', 'Jan', 'Feb', 'Dec', 'Mar', 'Jan', 'Feb', 'Nov', 'Nov', 'Mar', 'Jan', 'Mar', 'Dec', 'Feb']

SSW_year_N96 = [1852, 1854, 1854, 1855, 1855, 1856, 1857, 1857, 1862, 1862, 1863, 1864, 1865, 1865, 1868, 1869, 1870, 1872, 1873, 1873, 1874, 1875, 1876, 1878, 1879, 1880, 1880, 1881, 1882, 1883, 1885, 1886, 1886, 1887, 1888, 1888, 1889, 1891, 1892, 1894, 1895, 1896, 1896, 1897, 1901, 1902, 1903, 1904, 1904, 1905, 1906, 1907, 1908, 1909, 1911, 1912, 1913, 1915, 1918, 1918, 1920, 1920, 1922, 1922, 1923, 1924, 1926, 1927, 1929, 1929, 1930, 1933, 1935, 1938, 1939, 1941, 1941, 1943, 1944, 1946, 1946, 1947, 1948, 1948, 1949, 1950, 1950, 1951, 1952, 1954, 1954, 1955, 1955, 1957, 1958, 1958, 1960, 1960, 1961, 1964, 1966, 1967, 1968, 1970, 1973, 1974, 1976, 1976, 1977, 1979, 1980, 1982, 1982, 1987, 1988, 1989, 1990, 1991, 1991, 1991, 1992, 1992, 1994, 1994, 1996, 1997, 2000, 2001, 2001, 2001, 2002, 2004, 2005, 2008, 2009, 2010, 2012, 2013, 2013, 2014, 2015, 2017, 2018, 2019, 2019, 2024, 2024, 2025, 2028, 2029, 2030, 2030, 2034, 2039, 2040, 2043, 2044, 2045, 2046, 2046, 2049, 2052, 2054, 2057, 2057, 2058, 2059, 2060, 2061, 2066, 2067, 2070, 2071, 2072, 2073, 2073, 2074, 2075, 2075, 2076, 2077, 2077, 2081, 2085, 2087, 2087, 2088, 2090, 2091, 2093, 2094, 2094, 2097, 2097, 2099, 2100, 2100, 2103, 2103, 2105, 2105, 2106, 2107, 2107, 2107, 2108, 2108, 2116, 2116, 2118, 2119, 2122, 2122, 2123, 2123, 2123, 2124, 2125, 2125, 2125, 2126, 2126, 2128, 2130, 2131, 2131, 2132, 2132, 2134, 2137, 2138, 2145, 2146, 2146, 2146, 2149, 2150, 2150, 2151, 2153, 2155, 2156, 2158, 2160, 2160, 2161, 2161, 2163, 2164, 2169, 2169, 2171, 2173, 2176, 2176, 2178, 2179, 2179, 2180, 2181, 2181, 2181, 2182, 2183, 2186, 2187, 2188, 2188, 2189, 2191, 2192, 2193, 2194, 2194, 2198, 2200, 2201, 2201, 2203, 2204, 2206, 2209, 2211, 2215, 2215, 2215, 2216, 2217, 2218, 2220, 2220, 2220, 2221, 2222, 2223, 2224, 2224, 2225, 2226, 2228, 2229, 2230, 2230, 2231, 2233, 2238, 2240, 2244, 2245, 2246, 2246, 2248, 2249, 2250, 2250, 2252, 2252, 2253, 2255, 2255, 2256, 2257, 2258, 2259, 2259, 2263, 2264, 2264, 2266, 2266, 2268, 2268, 2269, 2272, 2274, 2275, 2275, 2275, 2276, 2277, 2278, 2279, 2283, 2284, 2287, 2291, 2292, 2294, 2295, 2295, 2296, 2302, 2304, 2306, 2307, 2307, 2308, 2308, 2309, 2310, 2310, 2313, 2313, 2314, 2316, 2319, 2320, 2320, 2322, 2322, 2326, 2327, 2330, 2331, 2332, 2337, 2337, 2339, 2340, 2341, 2342, 2343, 2344, 2345, 2345, 2346, 2349]

SSW_N216 = ['Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Jan', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Dec', 'Feb', 'Jan', 'Feb', 'Feb', 'Nov', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan', 'Mar', 'Mar', 'Mar', 'Feb', 'Nov', 'Nov', 'Jan', 'Jan', 'Mar', 'Feb', 'Feb', 'Mar', 'Feb', 'Jan', 'Mar', 'Jan', 'Mar', 'Feb', 'Dec', 'Feb', 'Mar', 'Jan', 'Nov', 'Feb', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Mar', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Nov', 'Mar', 'Mar', 'Jan', 'Mar', 'Feb', 'Feb', 'Jan', 'Jan', 'Feb', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Nov', 'Feb', 'Jan', 'Jan', 'Feb', 'Dec', 'Feb', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Mar', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Jan', 'Jan', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Dec', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Nov', 'Jan', 'Feb', 'Feb', 'Mar', 'Dec', 'Jan', 'Mar', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Feb', 'Jan', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Mar', 'Feb']

SSW_year_N216 = [1852, 1855, 1857, 1858, 1860, 1860, 1861, 1862, 1862, 1863, 1863, 1864, 1865, 1866, 1870, 1872, 1875, 1880, 1884, 1885, 1889, 1891, 1891, 1892, 1893, 1893, 1897, 1898, 1898, 1899, 1899, 1900, 1902, 1909, 1911, 1915, 1916, 1918, 1919, 1922, 1922, 1923, 1924, 1926, 1927, 1928, 1929, 1929, 1930, 1934, 1934, 1935, 1936, 1937, 1937, 1939, 1939, 1942, 1942, 1943, 1944, 1944, 1946, 1947, 1948, 1951, 1952, 1953, 1954, 1957, 1963, 1965, 1966, 1967, 1967, 1967, 1969, 1970, 1970, 1971, 1973, 1974, 1977, 1980, 1980, 1981, 1982, 1985, 1988, 1989, 1990, 1992, 1993, 1995, 1998, 2002, 2003, 2010, 2012, 2013, 2015, 2017, 2017, 2020, 2020, 2020, 2021, 2022, 2023, 2023, 2023, 2025, 2029, 2029, 2030, 2032, 2035, 2038, 2039, 2040, 2041, 2042, 2048, 2048, 2048, 2049, 2050, 2051, 2053, 2053, 2054, 2056, 2058, 2059, 2059, 2060, 2061, 2064, 2070, 2071, 2072, 2073, 2075, 2077, 2078, 2079, 2081, 2081, 2082, 2083, 2084, 2084, 2085, 2087, 2087, 2089, 2090, 2091, 2093, 2094, 2096, 2099, 2099, 2101, 2102, 2103, 2104, 2105, 2107, 2108, 2109, 2111, 2113, 2114, 2115, 2118, 2121, 2122, 2123, 2123, 2124, 2127, 2127, 2129]

SSW_N96_counter = Counter(SSW_N96)
SSW_N216_counter = Counter(SSW_N216)



per_decadeN96, ErN96, running_N96, year_countN96 = DJF_warmings_per_decade(SSW_N96, SSW_year_N96, years_N96)
per_decadeN216, ErN216, running_N216, year_countN216 = DJF_warmings_per_decade(SSW_N216, SSW_year_N216, years_N216)

fig = plt.figure(figsize=(18, 9))
plt.plot(years_N96[0:-10], running_N96)
plt.xlim(years_N96[0], years_N96[-1])
#plt.ylim(-0.1,5)
#plt.yticks([0,1,2,3,4,5])
plt.title('DJF SSWs per year, GC3 pi control N96 1850-2349')
plt.ylabel('number of SSWs')
plt.xlabel('year')
plt.show()



ps = np.abs(np.fft.fft(running_N216))**2

time_step = 1
freqs = np.fft.fftfreq(running_N216.size, time_step)
keep = freqs>=0
ps = ps[keep]
freqs = freqs[keep]
idx = np.argsort(freqs)



fig = plt.figure(figsize=(18, 9)) 
gs = gridspec.GridSpec(2, 1)

plt1 = plt.subplot(gs[0,0])
plt.plot(years_N216[0:-10], running_N216)
plt.xlim(years_N216[0], years_N216[-1])
#plt.ylim(-0.1,5)
#plt.yticks([0,1,2,3,4,5])
plt.title('DJF SSWs per year, GC3 pi control N216')
plt.ylabel('number of SSWs')
plt.xlabel('year')


plt2 = plt.subplot(gs[1,0])
plt.title('SSW count per decade power spectrum')
plt.plot(freqs[idx][1:-1], ps[idx][1:-1])
plt.xlabel('frequency (1/years)')
plt.ylabel('power (numebr of SSWs)$^2$')
plt.show()
fig.savefig('../figures/SSWN216_power_spec.png')






X = np.arange(500)
Y = 3*np.sin(2*math.pi*X/50)



ps = np.abs(np.fft.fft(Y))**2

time_step = 1
freqs = np.fft.fftfreq(Y.size, time_step)
keep = freqs>=0
ps = ps[keep]
freqs = freqs[keep]
idx = np.argsort(freqs)


fig = plt.figure(figsize=(18, 9)) 
gs = gridspec.GridSpec(2, 1)

plt1 = plt.subplot(gs[0,0])
plt.plot(X, Y)
#plt.ylim(-0.1,5)
#plt.yticks([0,1,2,3,4,5])
plt.xlabel('year')


plt2 = plt.subplot(gs[1,0])
plt.title('power spectrum')
plt.plot(freqs[idx][1:-1], ps[idx][1:-1])
plt.xlabel('frequency (1/years)')
plt.ylabel('power')
plt.show()
fig.savefig('../figures/test_power_spec.png')






