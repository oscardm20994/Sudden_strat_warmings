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


months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']

def time_add(List, tname):
	for i in range(len(List)):
		List[i] = add_times(List[i], tname)
	return List
	

#read in UKESM ensemble members
fname = '/gws/nopw/j04/aopp/strat_clim_data/CMIP6/UKESM1/historical/zmzw/ZMZW*.nc'
cubes = iris.load(fname)#, callback=callback)
cubes = cubes.extract(iris.Constraint(air_pressure = 10))
cubes = time_add(cubes, 'time')

#reanalysis
fname = '/gws/nopw/j04/aopp/strat_clim_data/reanalysis/ERA-5/daily/u_10hpa_daily_era5_1979-2019.nc'
ERA5 = iris.load_cube(fname)#, callback=callback)
ERA5 = add_times(ERA5, 'time')

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


	mean_perdec = np.mean(year_count)
	ER_perdec = np.std(year_count)/np.sqrt(len(years))
	return mean_perdec, ER_perdec
	


def get_SSW_rate(cube, tname):
	#start and end years for data
	start_year = cube.coord('year').points[0]
	end_year = cube.coord('year').points[-1]
	years = np.arange(start_year, end_year, 1)
	
	#take zonal mean
	cube = cube.collapsed('longitude', iris.analysis.MEAN)
	
	#get 60N
	cube = cube.extract(iris.Constraint(latitude = 60))
	
	#count SSWs
	SSW, SSW_year, indexi, indexj = SSW_counter(cube, start_year, end_year, 0)

	per_season, Er_perseason = DJF_warmings_per_decade(SSW, SSW_year, years)

	return SSW, per_season, Er_perseason, SSW_year

def by_mon(SSW, SSW_year, years):
	months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']
	month_count = []
	Er = []
	for i in range(len(months)):
		month_count.append(mon_warmings_per_decade(SSW, SSW_year, years, months[i])[0])
		Er.append(mon_warmings_per_decade(SSW, SSW_year, years, months[i])[1])
	return month_count, Er





#count SSWs in ERA5, ERAinterim below
SSW_ERA, per_seasonERA, Er_per_seasonERA, SSW_yearERA = get_SSW_rate(ERA5, 'time')
ERA_mon, ERA_mon_Er = by_mon(SSW_ERA, SSW_yearERA, np.arange(1979,2019,1))


#ERA interim SSWs
SSW_ERA = ['Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Jan', 'Jan', 'Dec', 'Mar', 'Feb', 'Dec', 'Feb', 'Mar', 'Feb', 'Dec', 'Jan', 'Jan', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan']
SSW_year_ERA = [1979, 1980, 1981, 1981, 1984, 1985, 1987, 1987, 1988, 1989, 1998, 1999, 2000, 2001, 2001, 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2010, 2013]
SSW_ERA40 = ['Jan', 'Dec', 'Feb', 'Jan', 'Nov', 'Mar', 'Jan', 'Jan', 'Mar', 'Jan', 'Jan']
SSW_year_ERA40 = [1963, 1965, 1966, 1968, 1968, 1969, 1970, 1971, 1971, 1973, 1977]
SSW_ERA2 = ['Feb', 'Feb', 'Mar']
SSW_year_ERA2 = [2017, 2018, 2018]
SSW_year_ERA = [ int(x) for x in SSW_year_ERA ]
SSW_ERA = SSW_ERA40 + SSW_ERA + SSW_ERA2
SSW_year_ERA = SSW_year_ERA40 + SSW_year_ERA + SSW_year_ERA2
per_seasonERA, Er_per_seasonERA = DJF_warmings_per_decade(SSW_ERA, SSW_year_ERA, np.arange(1958,2019,1))
ERA_mon, ERA_mon_Er = by_mon(SSW_ERA, SSW_year_ERA, np.arange(1958,2019,1))




#count SSWs in each UKESM ensemble
per_UK = []
ErUK = []
all_months = []
all_years = []
for i in range(len(cubes)):
	Set = get_SSW_rate(cubes[i], 'time')
	per_UK.append(Set[1])
	ErUK.append(Set[2]) 
	all_months.append(Set[0])
	all_years.append(Set[3])

all_months = [['Mar', 'Jan', 'Dec', 'Feb', 'Feb', 'Jan', 'Nov', 'Nov', 'Feb', 'Dec', 'Feb', 'Jan', 'Nov', 'Feb', 'Jan', 'Nov', 'Mar', 'Nov', 'Nov', 'Feb', 'Jan', 'Feb', 'Mar', 'Dec', 'Feb', 'Mar', 'Nov', 'Dec', 'Mar', 'Nov', 'Mar', 'Dec', 'Mar', 'Mar', 'Jan', 'Nov', 'Mar', 'Mar', 'Mar', 'Nov', 'Jan', 'Feb', 'Feb', 'Nov', 'Feb', 'Feb', 'Mar', 'Jan', 'Jan', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Dec', 'Feb', 'Dec', 'Feb', 'Mar', 'Dec', 'Feb', 'Nov', 'Nov', 'Mar', 'Feb', 'Dec', 'Mar', 'Mar', 'Feb', 'Jan', 'Mar', 'Nov', 'Nov', 'Feb', 'Jan', 'Feb', 'Mar', 'Nov', 'Jan', 'Mar', 'Nov', 'Nov', 'Nov', 'Nov', 'Feb', 'Mar', 'Mar', 'Dec', 'Jan', 'Jan', 'Mar', 'Feb'], ['Mar', 'Dec', 'Feb', 'Jan', 'Nov', 'Jan', 'Feb', 'Jan', 'Mar', 'Feb', 'Nov', 'Jan', 'Nov', 'Nov', 'Jan', 'Dec', 'Mar', 'Feb', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Jan', 'Feb', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Nov', 'Mar', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Dec', 'Jan', 'Jan', 'Nov', 'Mar', 'Nov', 'Feb', 'Feb', 'Feb', 'Nov', 'Nov', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Jan', 'Feb', 'Feb', 'Dec', 'Mar', 'Jan', 'Feb', 'Mar', 'Mar', 'Mar', 'Mar', 'Dec', 'Jan', 'Nov', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Feb', 'Nov', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Dec', 'Feb', 'Mar', 'Dec', 'Jan', 'Mar', 'Jan', 'Feb', 'Mar', 'Feb', 'Feb'], ['Jan', 'Feb', 'Dec', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Jan', 'Mar', 'Mar', 'Feb', 'Dec', 'Jan', 'Nov', 'Mar', 'Dec', 'Nov', 'Feb', 'Jan', 'Mar', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Mar', 'Dec', 'Mar', 'Dec', 'Feb', 'Mar', 'Jan', 'Jan', 'Feb', 'Nov', 'Nov', 'Nov', 'Dec', 'Feb', 'Jan', 'Dec', 'Feb', 'Nov', 'Mar', 'Jan', 'Mar', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Mar', 'Nov', 'Dec', 'Mar', 'Nov', 'Feb', 'Mar', 'Feb', 'Feb', 'Mar', 'Mar', 'Dec', 'Mar', 'Mar', 'Nov', 'Jan', 'Feb', 'Dec', 'Feb', 'Dec', 'Mar', 'Feb', 'Feb', 'Mar', 'Nov', 'Jan', 'Jan', 'Jan', 'Mar', 'Feb', 'Jan', 'Feb', 'Jan'], ['Mar', 'Jan', 'Jan', 'Nov', 'Jan', 'Nov', 'Feb', 'Mar', 'Feb', 'Nov', 'Mar', 'Feb', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Nov', 'Feb', 'Mar', 'Jan', 'Feb', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Jan', 'Mar', 'Jan', 'Dec', 'Feb', 'Mar', 'Mar', 'Mar', 'Jan', 'Dec', 'Dec', 'Feb', 'Jan', 'Nov', 'Jan', 'Mar', 'Feb', 'Mar', 'Jan', 'Nov', 'Feb', 'Jan', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Jan', 'Jan', 'Feb', 'Mar', 'Feb', 'Jan', 'Feb', 'Nov', 'Jan', 'Mar', 'Nov', 'Dec', 'Feb', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Dec', 'Dec', 'Feb', 'Mar', 'Jan', 'Nov', 'Jan', 'Jan', 'Dec', 'Feb', 'Nov', 'Mar', 'Mar', 'Jan', 'Feb', 'Feb', 'Jan', 'Dec', 'Jan', 'Feb', 'Dec', 'Jan', 'Nov', 'Jan', 'Feb', 'Feb', 'Dec', 'Jan', 'Feb', 'Nov', 'Dec', 'Mar', 'Feb', 'Feb'], ['Feb', 'Mar', 'Nov', 'Jan', 'Jan', 'Nov', 'Mar', 'Nov', 'Jan', 'Nov', 'Jan', 'Mar', 'Nov', 'Dec', 'Mar', 'Jan', 'Feb', 'Jan', 'Feb', 'Feb', 'Jan', 'Jan', 'Dec', 'Feb', 'Jan', 'Dec', 'Nov', 'Nov', 'Nov', 'Feb', 'Nov', 'Jan', 'Mar', 'Jan', 'Mar', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Nov', 'Nov', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Nov', 'Nov', 'Mar', 'Feb', 'Jan', 'Mar', 'Jan', 'Jan', 'Feb', 'Jan', 'Jan', 'Mar', 'Feb', 'Mar', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Jan', 'Feb', 'Nov', 'Nov', 'Nov', 'Mar', 'Mar', 'Nov', 'Jan', 'Jan', 'Feb', 'Nov', 'Nov', 'Jan', 'Mar', 'Jan', 'Nov', 'Nov', 'Jan', 'Mar', 'Dec', 'Nov', 'Mar', 'Nov', 'Dec', 'Mar', 'Dec', 'Mar', 'Nov', 'Feb', 'Feb'], ['Feb', 'Jan', 'Nov', 'Nov', 'Nov', 'Mar', 'Nov', 'Mar', 'Nov', 'Jan', 'Feb', 'Nov', 'Dec', 'Jan', 'Feb', 'Nov', 'Dec', 'Dec', 'Dec', 'Dec', 'Feb', 'Nov', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Dec', 'Feb', 'Nov', 'Mar', 'Dec', 'Feb', 'Feb', 'Nov', 'Mar', 'Jan', 'Mar', 'Mar', 'Nov', 'Mar', 'Nov', 'Feb', 'Feb', 'Mar', 'Nov', 'Jan', 'Jan', 'Mar', 'Nov', 'Jan', 'Feb', 'Mar', 'Feb', 'Mar', 'Mar', 'Dec', 'Mar', 'Feb', 'Jan', 'Dec', 'Mar', 'Mar', 'Mar', 'Mar', 'Feb', 'Dec', 'Feb', 'Jan', 'Feb', 'Nov', 'Jan', 'Feb', 'Mar', 'Nov', 'Mar', 'Jan', 'Feb', 'Mar', 'Mar', 'Mar', 'Mar', 'Nov', 'Feb', 'Feb', 'Nov', 'Jan', 'Feb', 'Mar', 'Mar'], ['Jan', 'Nov', 'Mar', 'Feb', 'Mar', 'Nov', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Jan', 'Mar', 'Dec', 'Feb', 'Mar', 'Dec', 'Nov', 'Feb', 'Jan', 'Feb', 'Jan', 'Jan', 'Feb', 'Feb', 'Feb', 'Dec', 'Mar', 'Nov', 'Jan', 'Nov', 'Feb', 'Jan', 'Jan', 'Jan', 'Feb', 'Dec', 'Feb', 'Mar', 'Feb', 'Mar', 'Feb', 'Jan', 'Jan', 'Mar', 'Nov', 'Mar', 'Nov', 'Feb', 'Jan', 'Feb', 'Nov', 'Jan', 'Dec', 'Mar', 'Nov', 'Mar', 'Dec', 'Feb', 'Feb', 'Dec', 'Mar', 'Jan', 'Feb', 'Mar', 'Dec', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Feb', 'Nov', 'Feb', 'Mar', 'Feb', 'Nov', 'Jan', 'Feb', 'Nov', 'Jan', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Dec', 'Feb', 'Jan', 'Feb', 'Nov', 'Nov', 'Jan', 'Feb', 'Nov', 'Jan', 'Mar', 'Nov', 'Feb', 'Dec', 'Dec', 'Feb', 'Mar', 'Dec'], ['Nov', 'Mar', 'Mar', 'Nov', 'Feb', 'Feb', 'Mar', 'Jan', 'Nov', 'Mar', 'Feb', 'Nov', 'Feb', 'Jan', 'Feb', 'Nov', 'Jan', 'Nov', 'Feb', 'Feb', 'Dec', 'Dec', 'Feb', 'Mar', 'Mar', 'Nov', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Dec', 'Mar', 'Dec', 'Dec', 'Nov', 'Jan', 'Feb', 'Feb', 'Feb', 'Mar', 'Nov', 'Mar', 'Feb', 'Jan', 'Jan', 'Mar', 'Dec', 'Feb', 'Dec', 'Jan', 'Mar', 'Jan', 'Feb', 'Mar', 'Mar', 'Mar', 'Nov', 'Nov', 'Dec', 'Dec', 'Nov', 'Feb', 'Nov', 'Jan', 'Feb', 'Feb', 'Mar', 'Feb', 'Nov', 'Mar', 'Feb', 'Jan', 'Nov', 'Feb', 'Nov', 'Feb', 'Dec', 'Feb', 'Jan', 'Feb', 'Feb', 'Feb', 'Mar', 'Mar', 'Nov', 'Mar', 'Nov', 'Jan', 'Mar', 'Jan', 'Mar', 'Nov', 'Nov', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Dec', 'Feb'], ['Jan', 'Feb', 'Feb', 'Feb', 'Mar', 'Nov', 'Jan', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Dec', 'Mar', 'Dec', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Feb', 'Dec', 'Dec', 'Feb', 'Nov', 'Feb', 'Jan', 'Feb', 'Nov', 'Feb', 'Feb', 'Feb', 'Nov', 'Dec', 'Mar', 'Nov', 'Mar', 'Feb', 'Jan', 'Mar', 'Feb', 'Mar', 'Nov', 'Dec', 'Feb', 'Nov', 'Jan', 'Feb', 'Mar', 'Nov', 'Feb', 'Mar', 'Nov', 'Feb', 'Feb', 'Feb', 'Dec', 'Jan', 'Nov', 'Nov', 'Jan', 'Dec', 'Nov', 'Feb', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan', 'Jan', 'Jan', 'Nov', 'Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Feb', 'Feb', 'Dec', 'Dec', 'Nov', 'Nov', 'Nov', 'Jan', 'Feb', 'Jan', 'Dec', 'Feb', 'Mar', 'Jan', 'Feb', 'Dec', 'Feb', 'Jan', 'Jan', 'Dec', 'Mar', 'Nov']]
all_years = [[1853, 1854, 1857, 1858, 1860, 1862, 1862, 1864, 1866, 1869, 1871, 1873, 1873, 1874, 1876, 1879, 1880, 1880, 1884, 1885, 1887, 1889, 1893, 1896, 1902, 1902, 1902, 1904, 1905, 1905, 1907, 1908, 1908, 1911, 1915, 1916, 1917, 1918, 1919, 1919, 1920, 1921, 1922, 1922, 1923, 1924, 1926, 1928, 1929, 1929, 1930, 1931, 1934, 1934, 1937, 1939, 1940, 1940, 1940, 1947, 1947, 1949, 1952, 1956, 1958, 1960, 1960, 1965, 1966, 1969, 1969, 1969, 1971, 1972, 1973, 1973, 1974, 1975, 1982, 1984, 1985, 1987, 1991, 1992, 1995, 2004, 2005, 2006, 2007, 2009, 2011, 2012], [1857, 1858, 1858, 1863, 1863, 1866, 1867, 1870, 1870, 1876, 1876, 1878, 1878, 1880, 1883, 1884, 1885, 1886, 1887, 1887, 1890, 1892, 1893, 1893, 1895, 1896, 1897, 1897, 1902, 1902, 1906, 1908, 1911, 1914, 1915, 1916, 1917, 1917, 1919, 1920, 1924, 1925, 1926, 1926, 1928, 1929, 1931, 1932, 1933, 1935, 1937, 1937, 1938, 1938, 1940, 1942, 1945, 1951, 1956, 1957, 1958, 1960, 1962, 1964, 1964, 1964, 1967, 1968, 1970, 1971, 1973, 1973, 1974, 1975, 1976, 1976, 1976, 1980, 1981, 1981, 1987, 1988, 1988, 1989, 1990, 1992, 1996, 1996, 1996, 1997, 2001, 2001, 2003, 2003, 2006, 2007, 2007, 2009, 2011, 2013], [1850, 1850, 1852, 1853, 1854, 1857, 1859, 1860, 1863, 1866, 1869, 1870, 1874, 1876, 1877, 1880, 1881, 1882, 1884, 1887, 1890, 1890, 1891, 1893, 1894, 1895, 1897, 1897, 1898, 1900, 1904, 1905, 1907, 1907, 1909, 1911, 1912, 1913, 1914, 1917, 1918, 1918, 1920, 1922, 1922, 1922, 1923, 1926, 1927, 1929, 1930, 1932, 1934, 1935, 1938, 1939, 1945, 1945, 1948, 1949, 1950, 1955, 1957, 1959, 1959, 1959, 1960, 1963, 1966, 1968, 1968, 1969, 1970, 1972, 1974, 1976, 1979, 1979, 1980, 1983, 1985, 1986, 1988, 1989, 1991, 1992, 2000, 2001, 2004, 2004, 2006, 2008, 2010, 2012], [1853, 1854, 1855, 1855, 1856, 1857, 1859, 1862, 1867, 1867, 1868, 1869, 1869, 1872, 1872, 1873, 1876, 1881, 1883, 1887, 1888, 1890, 1890, 1893, 1894, 1894, 1894, 1895, 1895, 1898, 1899, 1899, 1902, 1903, 1905, 1907, 1909, 1910, 1911, 1912, 1914, 1914, 1917, 1918, 1923, 1925, 1928, 1929, 1930, 1932, 1933, 1935, 1935, 1937, 1939, 1939, 1941, 1943, 1944, 1944, 1944, 1945, 1946, 1946, 1946, 1948, 1949, 1949, 1951, 1952, 1952, 1955, 1955, 1955, 1956, 1957, 1957, 1958, 1959, 1959, 1960, 1961, 1963, 1965, 1966, 1967, 1969, 1969, 1971, 1974, 1977, 1980, 1980, 1981, 1982, 1984, 1985, 1989, 1990, 1991, 1994, 1997, 1998, 1999, 2001, 2003, 2005, 2007, 2008, 2008, 2008, 2008, 2010, 2010, 2011, 2012], [1850, 1852, 1853, 1854, 1858, 1858, 1860, 1860, 1861, 1862, 1864, 1865, 1868, 1870, 1870, 1872, 1873, 1874, 1874, 1876, 1879, 1880, 1883, 1888, 1890, 1891, 1891, 1896, 1897, 1898, 1901, 1904, 1905, 1909, 1912, 1915, 1916, 1917, 1918, 1920, 1920, 1927, 1930, 1931, 1937, 1938, 1939, 1939, 1941, 1944, 1945, 1947, 1948, 1949, 1950, 1952, 1953, 1954, 1956, 1956, 1959, 1961, 1963, 1965, 1968, 1972, 1973, 1976, 1977, 1977, 1978, 1979, 1980, 1982, 1983, 1984, 1985, 1985, 1987, 1989, 1991, 1992, 1995, 1995, 1997, 1999, 1999, 2000, 2000, 2001, 2001, 2003, 2003, 2006, 2008, 2009, 2011, 2013], [1851, 1853, 1853, 1854, 1858, 1859, 1859, 1860, 1862, 1863, 1864, 1864, 1865, 1867, 1867, 1867, 1868, 1873, 1875, 1877, 1879, 1879, 1880, 1882, 1882, 1884, 1887, 1890, 1892, 1894, 1898, 1905, 1909, 1913, 1915, 1919, 1919, 1924, 1927, 1930, 1932, 1934, 1936, 1936, 1938, 1938, 1940, 1943, 1944, 1946, 1946, 1947, 1948, 1948, 1948, 1949, 1949, 1953, 1954, 1954, 1955, 1956, 1957, 1957, 1960, 1961, 1963, 1966, 1967, 1969, 1970, 1970, 1975, 1976, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1984, 1984, 1989, 1996, 2000, 2002, 2010, 2010, 2012], [1850, 1851, 1852, 1853, 1854, 1861, 1863, 1864, 1864, 1867, 1868, 1874, 1874, 1875, 1875, 1876, 1877, 1877, 1878, 1880, 1880, 1883, 1885, 1888, 1889, 1890, 1892, 1892, 1892, 1894, 1896, 1898, 1899, 1900, 1901, 1901, 1905, 1906, 1912, 1914, 1914, 1915, 1917, 1918, 1918, 1919, 1921, 1922, 1926, 1927, 1930, 1931, 1935, 1936, 1936, 1936, 1938, 1942, 1943, 1944, 1950, 1950, 1952, 1953, 1958, 1959, 1960, 1962, 1962, 1963, 1964, 1966, 1966, 1969, 1969, 1971, 1972, 1976, 1977, 1980, 1982, 1986, 1986, 1987, 1987, 1988, 1990, 1990, 1991, 1991, 1993, 1994, 1997, 1997, 1997, 1999, 2001, 2001, 2002, 2003, 2008, 2008, 2008, 2013], [1854, 1855, 1856, 1856, 1857, 1859, 1859, 1860, 1860, 1861, 1864, 1864, 1865, 1869, 1869, 1869, 1870, 1870, 1872, 1875, 1876, 1878, 1878, 1878, 1880, 1882, 1883, 1884, 1886, 1887, 1888, 1889, 1896, 1898, 1901, 1906, 1908, 1908, 1909, 1914, 1915, 1915, 1916, 1919, 1922, 1924, 1924, 1927, 1927, 1928, 1932, 1934, 1940, 1940, 1941, 1942, 1947, 1949, 1950, 1951, 1952, 1952, 1954, 1956, 1958, 1958, 1961, 1964, 1966, 1966, 1967, 1969, 1972, 1972, 1973, 1973, 1974, 1976, 1977, 1978, 1980, 1981, 1984, 1984, 1985, 1985, 1989, 1990, 1991, 1991, 1994, 1995, 1995, 1998, 2001, 2002, 2004, 2007, 2009, 2011, 2012, 2013], [1854, 1857, 1860, 1861, 1864, 1864, 1865, 1865, 1867, 1868, 1871, 1871, 1873, 1873, 1875, 1878, 1879, 1879, 1882, 1885, 1887, 1888, 1889, 1893, 1894, 1896, 1898, 1901, 1901, 1902, 1904, 1906, 1908, 1910, 1910, 1911, 1912, 1913, 1915, 1917, 1918, 1921, 1922, 1924, 1924, 1927, 1929, 1931, 1931, 1931, 1932, 1933, 1934, 1935, 1939, 1941, 1943, 1948, 1948, 1949, 1951, 1953, 1953, 1954, 1955, 1957, 1958, 1962, 1963, 1964, 1967, 1967, 1967, 1970, 1973, 1973, 1975, 1975, 1977, 1978, 1979, 1980, 1980, 1981, 1986, 1987, 1991, 1993, 1996, 1998, 1998, 2001, 2001, 2003, 2004, 2007, 2009, 2010, 2011, 2011]]


#analyse SSWs seperated by month
ensemble_mon = np.zeros([5,9])
ensemble_mon_Er = np.zeros([5,9])
for i in range(9):
	month_count, Er = by_mon(all_months[i], all_years[i], np.arange(1850,2015,1))
	ensemble_mon[:,i] = month_count
	ensemble_mon_Er[:,i] = Er
ensemble_mon = np.mean(ensemble_mon, axis = 1)
ensemble_mon_Er = np.sqrt(np.sum(ensemble_mon_Er**2, axis = 1))*ensemble_mon


#take mean over the ensemble
ensemble_mean_perdecade = np.mean(np.array(per_UK))
ensemble_er = np.sqrt(np.sum(np.array(ErUK)**2))*ensemble_mean_perdecade







#count SSWs in PI control run
SSW_UKESM = ['Jan', 'Dec', 'Feb', 'Mar', 'Dec', 'Dec', 'Jan', 'Dec', 'Nov', 'Feb', 'Mar', 'Dec', 'Mar', 'Mar', 'Dec', 'Feb', 'Nov', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Jan', 'Nov', 'Feb', 'Nov', 'Mar', 'Mar', 'Dec', 'Feb', 'Feb', 'Jan', 'Nov', 'Jan', 'Nov', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Nov', 'Feb', 'Feb', 'Mar', 'Mar', 'Mar', 'Mar', 'Feb', 'Jan', 'Nov', 'Feb', 'Feb', 'Nov', 'Mar', 'Feb', 'Jan', 'Jan', 'Dec', 'Jan', 'Mar', 'Jan', 'Mar', 'Nov', 'Feb', 'Feb', 'Mar', 'Feb', 'Feb', 'Jan', 'Feb', 'Nov', 'Mar', 'Mar', 'Nov', 'Mar', 'Mar', 'Jan', 'Feb', 'Dec', 'Jan', 'Nov', 'Nov', 'Jan', 'Feb', 'Nov', 'Feb', 'Feb', 'Mar', 'Dec', 'Mar', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Mar', 'Feb', 'Mar', 'Feb', 'Feb', 'Feb', 'Nov', 'Jan', 'Jan', 'Feb', 'Feb', 'Feb', 'Nov', 'Feb', 'Mar', 'Nov', 'Feb', 'Mar', 'Jan', 'Feb', 'Nov', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Nov', 'Dec', 'Dec', 'Mar', 'Nov', 'Feb', 'Dec', 'Jan', 'Feb', 'Nov', 'Jan', 'Feb', 'Feb', 'Mar', 'Mar', 'Jan', 'Mar', 'Dec', 'Feb', 'Dec', 'Mar', 'Feb', 'Jan', 'Feb', 'Nov', 'Feb', 'Mar', 'Mar', 'Nov', 'Mar', 'Feb', 'Feb', 'Mar', 'Jan', 'Mar', 'Nov', 'Jan', 'Mar', 'Mar', 'Jan', 'Feb', 'Mar', 'Jan', 'Dec', 'Feb', 'Mar', 'Feb', 'Jan', 'Jan', 'Nov', 'Jan', 'Mar', 'Jan', 'Dec', 'Mar', 'Feb', 'Mar', 'Feb', 'Nov', 'Jan', 'Jan', 'Mar', 'Feb', 'Jan', 'Nov', 'Feb', 'Jan', 'Nov', 'Nov', 'Nov', 'Feb', 'Mar', 'Feb', 'Mar', 'Mar', 'Mar', 'Nov', 'Jan', 'Jan', 'Nov', 'Mar', 'Jan', 'Feb', 'Nov', 'Mar']
SSW_year_UKESM = [2061, 2063, 2063, 2063, 2064, 2067, 2067, 2069, 2072, 2076, 2079, 2082, 2083, 2085, 2089, 2092, 2093, 2094, 2095, 2096, 2100, 2101, 2107, 2110, 2111, 2111, 2113, 2114, 2115, 2116, 2118, 2121, 2121, 2122, 2122, 2123, 2126, 2129, 2130, 2131, 2131, 2132, 2135, 2135, 2140, 2142, 2143, 2144, 2145, 2146, 2149, 2150, 2151, 2152, 2153, 2154, 2156, 2159, 2162, 2164, 2165, 2169, 2169, 2170, 2171, 2173, 2174, 2175, 2179, 2179, 2179, 2180, 2182, 2182, 2183, 2187, 2191, 2192, 2194, 2197, 2201, 2205, 2209, 2209, 2212, 2213, 2214, 2215, 2217, 2218, 2219, 2220, 2222, 2223, 2227, 2228, 2228, 2231, 2234, 2238, 2245, 2246, 2246, 2252, 2254, 2257, 2258, 2259, 2259, 2260, 2260, 2265, 2266, 2267, 2269, 2270, 2270, 2272, 2277, 2277, 2280, 2280, 2280, 2282, 2283, 2285, 2285, 2286, 2288, 2289, 2291, 2293, 2294, 2294, 2297, 2297, 2300, 2312, 2313, 2315, 2315, 2316, 2321, 2322, 2325, 2329, 2329, 2330, 2332, 2335, 2335, 2337, 2339, 2341, 2341, 2343, 2346, 2346, 2347, 2349, 2352, 2354, 2356, 2356, 2360, 2362, 2365, 2368, 2369, 2374, 2376, 2376, 2379, 2379, 2380, 2381, 2392, 2393, 2393, 2395, 2395, 2398, 2404, 2404, 2405, 2415, 2416, 2418, 2421, 2422, 2425, 2425, 2426, 2428, 2429, 2430, 2439, 2441, 2443, 2444, 2446, 2447, 2448, 2450, 2450, 2456, 2458] 
control_mean, control_Er = DJF_warmings_per_decade(SSW_UKESM, SSW_year_UKESM, np.arange(2060,2461,1))
control_Er = control_Er/np.sqrt(400)




index = np.arange(2)
fig = plt.subplots(figsize=(8,6), dpi = 400)
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3]) 
ax = plt.subplot(gs[0])
bar_width = 0.5
#warm1 = ax.bar(index[0:-1], np.sort(per_decadeUK), bar_width, yerr=ErUK,  capsize = 10, ecolor = 'black', color = 'r', linewidth = 2.0, label = 'Historical')
warm1 = ax.bar(1, per_seasonERA, bar_width, yerr=Er_per_seasonERA,  capsize = 7, ecolor = 'black', color = 'orange', linewidth = 2.0, label = 'ERA-interim')
warm2 = ax.bar(0, ensemble_mean_perdecade, bar_width, yerr=ensemble_er, capsize = 7, ecolor = 'black', color = 'blue', linewidth = 2.0, label = 'UKESM')
plt.ylabel('SSWs per season', fontsize = 13)
plt.xticks([0,1],[' ',' '])
plt.yticks(np.arange(0,0.9,0.1), fontsize = 13)
plt.xlim(-0.5,1.5)
#plt.ylim(0,0.6) 

ax1 = plt.subplot(gs[1])
index = np.arange(5)
bar_width = 0.25
warm1 = ax1.bar(index, ensemble_mon, bar_width, yerr=ensemble_mon_Er,  capsize = 7, ecolor = 'black', color = 'blue', linewidth = 2.0)
warm2 = ax1.bar(index + bar_width, ERA_mon, bar_width, yerr=ERA_mon_Er, capsize = 7, ecolor = 'black', color = 'orange', linewidth = 2.0)
ax1.set_xticks(index + bar_width/2)
ax1.set_xticklabels((months), fontsize = 13)
ax1.legend((warm1[0],warm2[0]),('UKESM1','ERA-Interim'), fontsize = 13, loc = 2)
plt.yticks(np.arange(6), fontsize = 13)
plt.xlim(-0.25,4.5)
plt.ylim(0,0.3) 
plt.yticks(np.arange(0,0.35,0.05), fontsize = 13)
plt.tight_layout()
plt.show()

fig[0].savefig('../../figures/Historical_SSWs_UKESM_ERAinterim_for_ACSIS.png', dpi = 200)










