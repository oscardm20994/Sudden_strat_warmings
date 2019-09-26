import matplotlib.pyplot as plt
import numpy as np
from GC3_diagnostics_functions import *
from scipy.optimize import curve_fit
import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import matplotlib.mlab as mlab
from scipy.stats import norm
from decimal import Decimal

#define count warmings per year
def warmings_per_decade(SSW, SSW_year, years):
	SSW_year2 = SSW_year[:]
	year_count = np.empty(0)
	running = np.empty(0)
	#### count number per decade
	SSW_year_counter = Counter(SSW_year2)
	for i in range(len(years)):
		year_count = np.append(year_count, SSW_year_counter[years[i]])
	print len(years)
	for i in range(len(years)-10):
		running = np.append(running, np.sum(year_count[i:i+10]))
	print running
	mean_perdec = np.mean(running)
	ER_perdec = np.std(running)
	#print running
	return mean_perdec, ER_perdec, running, year_count

def warmings_per_winter(year_count):
	per_season = np.mean(year_count)
	Er = np.std(year_count)/np.sqrt(len(year_count))
	
	return per_season, Er

#define count warmings per year from a given month
def mon_warmings_per_decade(SSW, SSW_year, years, mon):
	SSW_year2 = SSW_year[:]
	year_count = np.empty(0)
	running = np.empty(0)
	#### filter to warmings of a single month
	for i in range(len(SSW)):
		if SSW[i] != str(mon):
			SSW_year2[i] = 0
	SSW_year2 = filter(lambda a: a != 0, SSW_year2)
	#### count number per decade
	SSW_year_counter = Counter(SSW_year2)
	for i in range(len(years)):
		year_count = np.append(year_count, SSW_year_counter[years[i]])
	return np.mean(year_count), np.std(year_count)/np.sqrt(len(year_count))
	

#define the reanalysis period (58 years)
rean_per = len(years_ERA)
## define monte carlo algorithm that picks random reanalysis periods (58 years) from model runs and returns number of SSWs period
def monte_carlo(year_count, SSW_year, years):
	SSW_perseason_rean = np.empty(0)
	for i in range(len(years)-rean_per):
		print i
		start = i
		period = year_count[start:start+rean_per]
		SSW_perseason_rean = np.append(SSW_perseason_rean, np.mean(period))
	return SSW_perseason_rean
	

def monte_carlo_per_month(SSW, SSW_year, years, month):
	SSW_year2 = SSW_year[:]
	year_count = np.empty(0)
	running = np.empty(0)
	#### filter to warmings of a single month
	for i in range(len(SSW)):
		if SSW[i] != str(month):
			SSW_year2[i] = 0
	SSW_year2 = filter(lambda a: a != 0, SSW_year2)
	print SSW_year2
	year_count_month = warmings_per_decade(SSW, SSW_year2, years)[3]
	SSW_per_rean_month = monte_carlo(year_count_month, SSW_year2, years)
	
	return SSW_per_rean_month
		


###### code that takes a list of SSW years and months and perfomrs a monte carlo style algorithm to construct a probability distribution of number of
###### warmings per 50 year period which the reanalysis dataset length

## SSW arrays from UKESM, HadgemGC3.1 and reanalysis datasets
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



####proudcue synthetic SSW_yearcount data for each model
def synthetic_SSW(year_count):
	n_0 = np.count_nonzero(year_count == 0)
	n_1 = np.count_nonzero(year_count == 1)
	n_2 = np.count_nonzero(year_count == 2)
	n_3 = np.count_nonzero(year_count == 3)
	E = np.mean(year_count)
	p_0 = (E**0)*np.exp(-1*E)
	p_1 = (E**1)*np.exp(-1*E)
	p_2 = (E**2)*np.exp(-1*E)/2
	p_3 = (E**3)*np.exp(-1*E)/6

	tot = p_0 + p_1 + p_2 + p_3
	p_0 = p_0 + 1-tot
	year_count_synth = np.random.choice(np.array([0,1,2,3]), (len(year_count)), p = np.array([p_0, p_1, p_2, p_3]))
	
	return year_count_synth

year_count_synthN96 = synthetic_SSW(year_countN96)
year_count_synthN216 = synthetic_SSW(year_countN216)
year_count_synthUKESM = synthetic_SSW(year_countUKESM)








###combine reanalysis datests
SSW_ERA = SSW_ERA40 + SSW_ERA + SSW_ERA2
SSW_year_ERA = SSW_year_ERA40 + SSW_year_ERA + SSW_year_ERA2

#year arrays for different periods of run
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

start_yearERA = 1959 
end_yearERA = 2018
years_ERA = np.arange(start_yearERA, end_yearERA, 1)

### call function from script GC3_diagnostics_functions.py to count number of warmings for each year in each dataset (year_count)
per_decadeERA, ErERA, running_ERA, year_countERA = warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA)
per_decadeN96, ErN96, running_N96, year_countN96 = warmings_per_decade(SSW_N96, SSW_year_N96, years_N96)
per_decadeN216, ErN216, running_N216, year_countN216 = warmings_per_decade(SSW_N216, SSW_year_N216, years_N216)
per_decadeUKESM, ErUKESM, running_UKESM, year_countUKESM = warmings_per_decade(SSW_UKESM, SSW_year_UKESM, years_UKESM)

##### find mean number of SSWs per winter
N96_per_season, Er_per_seasonN96 = warmings_per_winter(year_countN96) 
N216_per_season, Er_per_seasonN216 = warmings_per_winter(year_countN216) 
UKESM_per_season, Er_per_seasonUKESM = warmings_per_winter(year_countUKESM) 
ERA_per_season, Er_per_seasonERA = warmings_per_winter(year_countERA) 

series = np.array([N96_per_season, N216_per_season, UKESM_per_season, ERA_per_season])
Er_series = np.array([Er_per_seasonN96, Er_per_seasonN216, Er_per_seasonUKESM, Er_per_seasonERA])
index = np.arange(4) + 0.5
fig, ax = plt.subplots(figsize=(18,12))
bar_width = 0.5
plt.bar(index, series, bar_width, yerr=Er_series,  capsize = 10, ecolor = 'black', color = ('r','g','b','cyan'), linewidth = 2.0)
ax.set_xticks(index)
ax.set_xticklabels((['N96','N216','UKESM','ERA']), fontsize = 20)
plt.title('Mean SSWs per winter', fontsize = 30)
plt.ylabel('SSWs per winter', fontsize = 20)
plt.yticks(np.arange(0,1,0.1), fontsize = 20)
plt.xlim(-0.25,4)
plt.ylim(0,1) 
plt.show()
fig.savefig('../../figures/SSWs_per_winter_bar.png')








#seperate pdfs by month
months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']

N96_set = np.empty([0,0])
UKESM_set = np.empty([0,0])
N216_set = np.empty([0,0])


SSW_per_reanN96_Jan = monte_carlo_per_month(SSW_N96, SSW_year_N96, years_N96, 'Jan') 





## plot composites for each month seperately for ERA and N96 data
#levels = np.arange(-3.2,np.max(weaks_mean_UKESM)+2,0.5)

binwidth = 0.2
ERA_count = sum(year_countERA)
fig, ax = plt.subplots(figsize=(18, 18))
gs = gridspec.GridSpec(len(months), 3)
for i in range(len(months)):
	SSW_per_reanN96 = monte_carlo_per_month(SSW_N96, SSW_year_N96, years_N96, months[i])
	SSW_per_reanN216 = monte_carlo_per_month(SSW_N216, SSW_year_N216, years_N216, months[i])
        SSW_per_reanUKESM = monte_carlo_per_month(SSW_UKESM, SSW_year_UKESM, years_UKESM, months[i])
	
	
	
### plot composites for all months
        plt1 = plt.subplot(gs[i,0])
	plt.hist(SSW_per_reanN96, bins=np.arange(min(SSW_per_reanN96), max(SSW_per_reanN96) + binwidth, binwidth), color = 'r', alpha = 0.5, density = True, label = 'HadgemGC3.1 N96')
	bins=np.arange(min(SSW_per_reanN96), max(SSW_per_reanN96) + binwidth, binwidth)
	(mu, sigma) = norm.fit(SSW_per_reanN96)
	y = mlab.normpdf( bins, mu, sigma)
	plt.plot(bins+0.5, y, color = 'r')
	if i == 0:
                plt.title('(Hadgem GC3.1 N96)'+ str(months[i]), fontsize = 15)
        else:
             	plt.title(str(months[i]), fontsize = 15)
		
	plt2 = plt.subplot(gs[i,1])
	plt.hist(SSW_per_reanUKESM, bins=np.arange(min(SSW_per_reanUKESM), max(SSW_per_reanUKESM) + binwidth, binwidth), color = 'g', alpha = 0.5, density = True, label = 'UKESM')
	bins=np.arange(min(SSW_per_reanUKESM), max(SSW_per_reanUKESM) + binwidth, binwidth)
	(mu, sigma) = norm.fit(SSW_per_reanUKESM)
	y = mlab.normpdf( bins, mu, sigma)
	plt.plot(bins+0.5, y, color = 'g')
	if i == 0:
                plt.title(' (UKESM)'+ str(months[i]), fontsize = 15)
        else:
             	plt.title(str(months[i]), fontsize = 15)
		
        plt3 = plt.subplot(gs[i,2])
	plt.hist(SSW_per_reanN216, bins=np.arange(min(SSW_per_reanN216), max(SSW_per_reanN216) + binwidth, binwidth), color = 'b', alpha = 0.5, density = True, label = 'HadgemGC3.1 N216')
	bins=np.arange(min(SSW_per_reanN216), max(SSW_per_reanN216) + binwidth, binwidth)
	(mu, sigma) = norm.fit(SSW_per_reanN216)
	y = mlab.normpdf( bins, mu, sigma)
	plt.plot(bins+0.5, y, color = 'b')
	if i == 0:
        	plt.title(' (Hadgem GC3.1 N216)'+ str(months[i]), fontsize = 15)
        else:
             	plt.title(str(months[i]), fontsize = 15)
		
plt.show()




































#call monte carlo function above for real and synthetic data
SSW_perseason_reanN96 = monte_carlo(year_countN96, SSW_year_N96, years_N96)
SSW_perseason_reanUKESM = monte_carlo(year_countUKESM, SSW_year_UKESM, years_UKESM)
SSW_perseason_reanN216 = monte_carlo(year_countN216, SSW_year_N216, years_N216)

SSW_perseason_reanN96_synth = monte_carlo(year_count_synthN96, SSW_year_N96, years_N96)
SSW_perseason_reanUKESM_synth = monte_carlo(year_count_synthUKESM, SSW_year_UKESM, years_UKESM)
SSW_perseason_reanN216_synth = monte_carlo(year_count_synthN216, SSW_year_N216, years_N216)



binwidth = 0.017
fig, ax = plt.subplots(figsize=(18, 18))
gs = gridspec.GridSpec(3, 2)

plt1 = plt.subplot(gs[0,0])
plt.hist(SSW_perseason_reanN96, bins=np.arange(min(SSW_perseason_reanN96), max(SSW_perseason_reanN96) + binwidth, binwidth), color = 'r', alpha = 0.4, density = True, label = 'HadgemGC3.1 N96')
bins=np.arange(min(SSW_perseason_reanN96), max(SSW_perseason_reanN96) + binwidth, binwidth)
(mu, sigma) = norm.fit(SSW_perseason_reanN96)
y = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, y, color = 'r')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black')
plt.xlabel('SSWs per season', fontsize = 25)
plt.xlim(0,1)
plt.legend(fontsize = 15)

plt2 = plt.subplot(gs[0,1])
plt.hist(SSW_perseason_reanN96_synth, bins=np.arange(min(SSW_perseason_reanN96_synth), max(SSW_perseason_reanN96_synth) + binwidth, binwidth), color = 'r', alpha = 0.4, density = True, label = 'HadgemGC3.1 N96 (random)')
bins=np.arange(min(SSW_perseason_reanN96_synth), max(SSW_perseason_reanN96_synth) + binwidth, binwidth)
(mu, sigma) = norm.fit(SSW_perseason_reanN96_synth)
y = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, y, color = 'r')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black')
plt.xlabel('SSWs per season', fontsize = 25)
plt.xlim(0,1)
plt.legend(fontsize = 15)

plt3 = plt.subplot(gs[1,0])
plt.hist(SSW_perseason_reanN216, bins=np.arange(min(SSW_perseason_reanN216), max(SSW_perseason_reanN216) + binwidth, binwidth), color = 'g', alpha = 0.4, density = True, label = 'HadgemGC3.1 N216')
bins=np.arange(min(SSW_perseason_reanN216), max(SSW_perseason_reanN216) + binwidth, binwidth)
(mu, sigma) = norm.fit(SSW_perseason_reanN216)
y = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, y, color = 'g')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black')
plt.xlabel('SSWs per season', fontsize = 25)
plt.xlim(0,1)
plt.legend(fontsize = 15)

plt4 = plt.subplot(gs[1,1])
plt.hist(SSW_perseason_reanN216_synth, bins=np.arange(min(SSW_perseason_reanN216_synth), max(SSW_perseason_reanN216_synth) + binwidth, binwidth), color = 'g', alpha = 0.4, density = True, label = 'HadgemGC3.1 N216 (random)')
bins=np.arange(min(SSW_perseason_reanN216_synth), max(SSW_perseason_reanN216_synth) + binwidth, binwidth)
(mu, sigma) = norm.fit(SSW_perseason_reanN216_synth)
y = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, y, color = 'g')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black')
plt.xlabel('SSWs per season', fontsize = 25)
plt.xlim(0,1)
plt.legend(fontsize = 15)

plt5 = plt.subplot(gs[2,0])
plt.hist(SSW_perseason_reanUKESM, bins=np.arange(min(SSW_perseason_reanUKESM), max(SSW_perseason_reanUKESM) + binwidth, binwidth), color = 'b', alpha = 0.4, density = True, label = 'UKESM1 N96')
bins=np.arange(min(SSW_perseason_reanUKESM), max(SSW_perseason_reanUKESM) + binwidth, binwidth)
(mu, sigma) = norm.fit(SSW_perseason_reanUKESM)
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, color = 'b')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black')
plt.xlabel('SSWs per season', fontsize = 25)
plt.xlim(0,1)
plt.legend(fontsize = 15)

plt6 = plt.subplot(gs[2,1])
plt.hist(SSW_perseason_reanUKESM_synth, bins=np.arange(min(SSW_perseason_reanUKESM_synth), max(SSW_perseason_reanUKESM_synth) + binwidth, binwidth), color = 'b', alpha = 0.4, density = True, label = 'UKESM (random)')
bins=np.arange(min(SSW_perseason_reanUKESM_synth), max(SSW_perseason_reanUKESM_synth) + binwidth, binwidth)
(mu, sigma) = norm.fit(SSW_perseason_reanUKESM_synth)
y = mlab.normpdf( bins, mu, sigma)
plt.plot(bins, y, color = 'b')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black')
plt.xlabel('SSWs per season', fontsize = 25)
plt.xlim(0,1)
plt.legend(fontsize = 15)

plt.show()

fig.savefig('../../figures/SSWs_perseason_RA_period')
















