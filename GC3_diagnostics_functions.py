import iris
import iris.coord_categorisation as coord_cat
import numpy as np
from collections import Counter
import sys
import matplotlib.pyplot as plt
import sys
#import vor_fast
#import vor_fast_setup
from matplotlib import gridspec
from mpl_toolkits.basemap import Basemap


###Functions called from 'GC3_diagnostics.py' and 'GC3_diagnostics_climatology.py'

#averages over longitude, month and season (to produce climatology contour plot)
def U_averaging(U):
    #average over longitudes
    U_lats = U.collapsed('longitude', iris.analysis.MEAN) 
    #average over months
    U_monthly = U_lats.aggregated_by('month', iris.analysis.MEAN)
    #average over seasons
    U_season = U_lats.aggregated_by('clim_season', iris.analysis.MEAN)
    
    return U_season, U_monthly

#constrain to 10hPa pressure level
def ten_hpa_U(U):
    U_lats = U.collapsed('longitude', iris.analysis.MEAN)
    constraint = iris.Constraint(Pressure=10.)#, latitude =6.00000038e+01)
    
    U_ten_hpa = U_lats.extract(constraint)
    #U_ten_hpa = U_lat_lon.extract(pressure_constraint)
    
    return U_ten_hpa

#constratin latitude to 60N    
def correct_latitude60(cell):
    """Return True for values that are in the correct range."""
    return    5.94444504e+01 < cell < 6.055556106e+01  
lat_constraint=iris.Constraint(latitude = correct_latitude60)


def correct_latitude60m(cell):
    """Return True for values that are in the correct range."""
    return      -6.055556106e+01 < cell < -5.94444504e+01 
lat_constraintm=iris.Constraint(latitude = correct_latitude60m)

def correct_latitudeNAON(cell):
    """Return True for values that are in the correct range."""
    return     55.0 < cell < 80.0
lat_constraintNAON=iris.Constraint(latitude = correct_latitudeNAON)

def correct_latitudeNAOS(cell):
    """Return True for values that are in the correct range."""
    return     20.0 < cell < 55.0
lat_constraintNAOS=iris.Constraint(latitude = correct_latitudeNAOS)

def correct_longitudeNAO1(cell):
    """Return True for values that are in the correct range."""
    return    0 < cell < 30.0 or 280.0 < cell
lon_constraintNAO1=iris.Constraint(longitude = correct_longitudeNAO1)

def correct_lat_EOF(cell):
	return 20 < cell < 80
lat_constraint_EOF = iris.Constraint(latitude = correct_lat_EOF)

def correct_lon_EOF(cell):
    return    260 < cell < 360 #or 0 < cell < 20
lon_constraint_EOF = iris.Constraint(longitude = correct_lon_EOF)


#extract annual seasonal mean time series

def seasonal_series(U, season, start_year, end_year):
	if season == 'djf':
		U_seasons = iris.cube.CubeList([])
		Ud = U.extract(iris.Constraint(month = 'Dec'))
		Ujf = U.extract(iris.Constraint(month = ['Jan', 'Feb']))
		for i in range(start_year, end_year):
			dummy2 = iris.cube.CubeList([Ud.extract(iris.Constraint(year = i)),Ujf.extract(iris.Constraint(year = i+1))])
			U_seasons.append(dummy2.concatenate_cube())
		U_season_series = np.empty(0)
		for i in range(len(U_seasons)):		
			U_seasons[i] = U_seasons[i].collapsed('t', iris.analysis.MEAN)
			U_season_series = np.append(U_season_series, U_seasons[i].data)


	else:
		U_season = U.extract(iris.Constraint(clim_season = season))
		U_season_series = U_season.aggregated_by('year', iris.analysis.MEAN)
	return U_season_series



#extract each month from years 
def monthly_series(U, month):
    U_month = U.extract(iris.Constraint(month = month))
    #U_month_series = U_month.aggregated_by('year', iris.analysis.MEAN)
    return U_month




def SSW_counter(U_strat, start_year, end_year, thresh):#, GPH):
#identify SSWs from U_strat cube. defition of an SSW: U becomes negative at 60N 10hPa in months NDJFM (20 days after an SSW, another can occur)

	#restrict months
	U_SSW_seasons, SSW, final, indexi, indexj, SSW_year = [], [], [], [], [], []
	U_strat_ND = U_strat.extract(iris.Constraint(month = ['Nov', 'Dec']))
	U_strat_JFMA = U_strat.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr']))


	#GPH_ND = GPH.extract(iris.Constraint(month = ['Nov', 'Dec']))
	#GPH_JFMA =GPH.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr']))
	

	U_SSW_seasons.append(U_strat_JFMA.extract(iris.Constraint(year = start_year)))
	for i in range(start_year, end_year-1):
		print i
		dummy2 = iris.cube.CubeList([U_strat_ND.extract(iris.Constraint(year = i)),U_strat_JFMA.extract(iris.Constraint(year = i+1))])
		U_SSW_seasons.append(dummy2.concatenate_cube())
		print U_SSW_seasons
		#dummy3 = iris.cube.CubeList([GPH_ND.extract(iris.Constraint(year = i)),GPH_JFMA.extract(iris.Constraint(year = i+1))])
		#GPH_seasons.append(dummy3.concatenate_cube())

	###plot anitmation of vortex
	#levels = np.arange(-40,80,5)
	#for i in range(len(GPH_seasons[33].data[:,0,0])):
	#	fig = plt.figure(figsize=(24,14))
	#	plt.title(str(i))
	#	m = Basemap(projection='spstere',lon_0=0,boundinglat=-45)
	#	lon,lat=np.meshgrid(GPH_seasons[33].coord('longitude').points,GPH_seasons[33].coord('latitude').points)
	#	x, y =m(lon,lat)
	#	m.drawcoastlines()
	#	m.contourf(x,y,GPH_seasons[33].data[i,:,:], origin='lower')
	#	plt.savefig('../figures/animations/'+str(i))
	#	plt.show()
	#sys.exit()
	
	
	for j in range(len(U_SSW_seasons)):
		print j	
		i=0
		while i < len(U_SSW_seasons[j].data):
			#print i,j
			if U_SSW_seasons[j].data[i] < thresh and U_SSW_seasons[j].data[i-1] >= thresh:
				SSW.append(str(U_SSW_seasons[j].coord('month').points[i]))
				SSW_year.append(U_SSW_seasons[j].coord('season_year').points[i])
				indexi.append(i)
				indexj.append(j)
                		#GPH_nh, lats_nh, xypoints = vor_fast_setup.setup(GPH_seasons[j].data[i,:,:] ,GPH_seasons[j].coord('latitude').points, GPH_seasons[j].coord('longitude').points,'NH')
                		#SSW_Kurt.append(vortex_calc(GPH_nh, lats_nh, GPH_seasons[j].coord('longitude').points, xypoints))
				final.append(i)
				for x in range(i, len(U_SSW_seasons[j].data)):
					if np.all(U_SSW_seasons[j].data[x:x+19] > thresh):
						i = x+19
						break

					elif x == len(U_SSW_seasons[j].data)-1:
						i = 1000000
						
					else:
						j=j	
			else:
				i = i+1
		

		
		
		#print 'SSWs before checking for warmings is', SSW 
		
		
		#Check for a fianl warming and remove from SSW list if found
		#f = np.where(U_SSW_seasons[j].data == final[-1])[0]
		#print f #U_SSW_seasons[j].data[f-1:f+1]
		
		if final != []:	
			for k in range(final[-1],len(U_SSW_seasons[j].data)-9):
				if np.all(U_SSW_seasons[j].data[k:k+9] > thresh):
					break
				elif k== len(U_SSW_seasons[j].data)-10: 
					#print 'final warming found'
					del SSW[-1]
					del SSW_year[-1]
                    			#del SSW_Kurt[-1]
					del indexi[-1]
					del indexj[-1]
					break
				else:
					k=k	
		else:	
			j = j
			
			
			
		### plot kurtosis 10 days either side of warming
	
		#print 'SSWs after checking for final warmings is', SSW
		#print U_SSW_seasons[j].data


		#fig = plt.figure(figsize=(18, 9))
		#plt.title('Winter westerly wind, 10hPa 60N and GPH vortex Kurtosis')
		#gs = gridspec.GridSpec(2, 1)
		
		#plt1 = plt.subplot(gs[0,0])
		#plt.plot(np.linspace(0,6,len(U_SSW_seasons[j].data)),U_SSW_seasons[j].data)
        	#plt.plot([0,6],[0,0])
		#plt.ylabel('U (ms -1)')
		#plt.xticks(np.arange(0,6,1),['01/11', '01/12', '01/01', '01/02', '01/03', '01/04'])
        	#plt.title(str(2011+j) + '/' + str(2012+j))
		

		#plt2 = plt.subplot(gs[1,0])
		#plt.plot(np.linspace(0,6,len(kurt)),kurt)
        	#plt.plot([0,6],[-0.1,-0.1])
		#plt.ylabel('Kurtosis')
		#plt.xticks(np.arange(0,6,1),['01/11', '01/12', '01/01', '01/02', '01/03', '01/04'])
		

		#fig.savefig('../figures/SSW_years/U_kurtosis_GC3_'+ str(2011+j) + '-' + str(2012+j) + '.png')
		
		
#plt.show()
			

	

	for i in range(len(SSW)):
		if SSW[i] == 'Apr':
			SSW_year[i] = 0
            		#SSW_Kurt = 0
			indexi[i] = 'xx'
			indexj[i] = 'xx'

	SSW_year = filter(lambda a: a != 0, SSW_year)
    	#SSW_Kurt= filter(lambda a: a != 0, SSW_Kurt)
    	SSW = filter(lambda a: a != 'Apr', SSW)
	indexi = filter(lambda a: a!= 'xx', indexi)
	indexj = filter(lambda a: a!= 'xx', indexj)
	#print SSW
	
	
	
	
	#kurt = []
	#for j in range(start_year, end_year):
	#	GPH_nh, lats_nh, xypoints = vor_fast_setup.setup(GPH_seasons[j].data,GPH_seasons[j].coord('latitude').points, GPH_seasons[j].coord('longitude').points,'NH')
			
	#	for p in range(len(GPH_nh[:,0,0])):
	#		kurt.append(vortex_calc(GPH_nh[p,:,:], lats_nh, GPH_seasons[j].coord('longitude').points, xypoints))
		
	print SSW
	#for x in index:
		#plt.plot(np.linspace(-10,10,len(kurt[x-20:x+20])), kurt[x-20:x+20])
		#plt.show()
		
	
	



	
	print SSW_year, len(SSW_year)
	#print SSW_index, len(SSW_index) 
	return SSW, SSW_year, indexi, indexj



def SSW_counter_SP(U_strat, start_year, end_year, thresh):#, GPH):
#identify SSWs from U_strat cube. defition of an SSW: U becomes negative at 60N 10hPa in months NDJFM (20 days after an SSW, another can occur)

	#restrict months
	U_SSW_seasons, SSW, final, indexi, indexj, SSW_year = [], [], [], [], [], []
	U_strat_ND = U_strat.extract(iris.Constraint(month = ['Jul','Aug','Sep','Oct','Nov',]))
	#U_strat_JFMA = U_strat.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', ]))


	#GPH_ND = GPH.extract(iris.Constraint(month = ['Nov', 'Dec']))
	#GPH_JFMA =GPH.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr']))
	

	U_SSW_seasons.append(U_strat_JFMA.extract(iris.Constraint(year = start_year)))
	for i in range(start_year, end_year-1):
		print i
		dummy2 = U_strat_ND.extract(iris.Constraint(year = i))
		U_SSW_seasons.append(dummy2)
		print U_SSW_seasons
		#dummy3 = iris.cube.CubeList([GPH_ND.extract(iris.Constraint(year = i)),GPH_JFMA.extract(iris.Constraint(year = i+1))])
		#GPH_seasons.append(dummy3.concatenate_cube())

	###plot anitmation of vortex
	#levels = np.arange(-40,80,5)
	#for i in range(len(GPH_seasons[33].data[:,0,0])):
	#	fig = plt.figure(figsize=(24,14))
	#	plt.title(str(i))
	#	m = Basemap(projection='spstere',lon_0=0,boundinglat=-45)
	#	lon,lat=np.meshgrid(GPH_seasons[33].coord('longitude').points,GPH_seasons[33].coord('latitude').points)
	#	x, y =m(lon,lat)
	#	m.drawcoastlines()
	#	m.contourf(x,y,GPH_seasons[33].data[i,:,:], origin='lower')
	#	plt.savefig('../figures/animations/'+str(i))
	#	plt.show()
	#sys.exit()
	
	
	for j in range(len(U_SSW_seasons)):
		print j	
		i=0
		while i < len(U_SSW_seasons[j].data):
			#print i,j
			if U_SSW_seasons[j].data[i] < thresh and U_SSW_seasons[j].data[i-1] >= thresh:
				SSW.append(str(U_SSW_seasons[j].coord('month').points[i]))
				SSW_year.append(U_SSW_seasons[j].coord('season_year').points[i])
				indexi.append(i)
				indexj.append(j)
                		#GPH_nh, lats_nh, xypoints = vor_fast_setup.setup(GPH_seasons[j].data[i,:,:] ,GPH_seasons[j].coord('latitude').points, GPH_seasons[j].coord('longitude').points,'NH')
                		#SSW_Kurt.append(vortex_calc(GPH_nh, lats_nh, GPH_seasons[j].coord('longitude').points, xypoints))
				final.append(i)
				for x in range(i, len(U_SSW_seasons[j].data)):
					if np.all(U_SSW_seasons[j].data[x:x+19] > thresh):
						i = x+19
						break

					elif x == len(U_SSW_seasons[j].data)-1:
						i = 1000000
						
					else:
						j=j	
			else:
				i = i+1
		

		
		
		#print 'SSWs before checking for warmings is', SSW 
		
		
		#Check for a fianl warming and remove from SSW list if found
		#f = np.where(U_SSW_seasons[j].data == final[-1])[0]
		#print f #U_SSW_seasons[j].data[f-1:f+1]
		
		if final != []:	
			for k in range(final[-1],len(U_SSW_seasons[j].data)-9):
				if np.all(U_SSW_seasons[j].data[k:k+9] > thresh):
					break
				elif k== len(U_SSW_seasons[j].data)-10: 
					#print 'final warming found'
					del SSW[-1]
					del SSW_year[-1]
                    			#del SSW_Kurt[-1]
					del indexi[-1]
					del indexj[-1]
					break
				else:
					k=k	
		else:	
			j = j
			
			
			
		### plot kurtosis 10 days either side of warming
	
		#print 'SSWs after checking for final warmings is', SSW
		#print U_SSW_seasons[j].data


		#fig = plt.figure(figsize=(18, 9))
		#plt.title('Winter westerly wind, 10hPa 60N and GPH vortex Kurtosis')
		#gs = gridspec.GridSpec(2, 1)
		
		#plt1 = plt.subplot(gs[0,0])
		#plt.plot(np.linspace(0,6,len(U_SSW_seasons[j].data)),U_SSW_seasons[j].data)
        	#plt.plot([0,6],[0,0])
		#plt.ylabel('U (ms -1)')
		#plt.xticks(np.arange(0,6,1),['01/11', '01/12', '01/01', '01/02', '01/03', '01/04'])
        	#plt.title(str(2011+j) + '/' + str(2012+j))
		

		#plt2 = plt.subplot(gs[1,0])
		#plt.plot(np.linspace(0,6,len(kurt)),kurt)
        	#plt.plot([0,6],[-0.1,-0.1])
		#plt.ylabel('Kurtosis')
		#plt.xticks(np.arange(0,6,1),['01/11', '01/12', '01/01', '01/02', '01/03', '01/04'])
		

		#fig.savefig('../figures/SSW_years/U_kurtosis_GC3_'+ str(2011+j) + '-' + str(2012+j) + '.png')
		
		
#plt.show()
			

	

	for i in range(len(SSW)):
		if SSW[i] == 'Apr':
			SSW_year[i] = 0
            		#SSW_Kurt = 0
			indexi[i] = 'xx'
			indexj[i] = 'xx'

	SSW_year = filter(lambda a: a != 0, SSW_year)
    	#SSW_Kurt= filter(lambda a: a != 0, SSW_Kurt)
    	SSW = filter(lambda a: a != 'Apr', SSW)
	indexi = filter(lambda a: a!= 'xx', indexi)
	indexj = filter(lambda a: a!= 'xx', indexj)
	#print SSW
	
	
	
	
	#kurt = []
	#for j in range(start_year, end_year):
	#	GPH_nh, lats_nh, xypoints = vor_fast_setup.setup(GPH_seasons[j].data,GPH_seasons[j].coord('latitude').points, GPH_seasons[j].coord('longitude').points,'NH')
			
	#	for p in range(len(GPH_nh[:,0,0])):
	#		kurt.append(vortex_calc(GPH_nh[p,:,:], lats_nh, GPH_seasons[j].coord('longitude').points, xypoints))
		
	print SSW
	#for x in index:
		#plt.plot(np.linspace(-10,10,len(kurt[x-20:x+20])), kurt[x-20:x+20])
		#plt.show()
		
	
	



	
	print SSW_year, len(SSW_year)
	#print SSW_index, len(SSW_index) 
	return SSW, SSW_year, indexi, indexj

















#add month, season and year metadata
def add_times(cube, time):
	coord_cat.add_month(cube, time, name='month')
	coord_cat.add_season(cube, time, name='clim_season')
	coord_cat.add_year(cube, time, name='year')
	coord_cat.add_day_of_year(cube, time, name='day_number')
	coord_cat.add_season_year(cube, time, name='season_year')

	return cube


#def SSW_in_years(U_strat, start_year, end_year):
#	SSW_which_year = SSW_counter(U_strat, start_year, end_year, 0)[1]
#	SSW_year_counter = Counter(SSW_years)
#	SSWs_years = SSW_year_counter[np.arange(start_year, end_year, 1)]
#	print SSW_year_counter
#	print SSWs_years
	
#	return SSW_year




#data for spaghetti plot
#constrain to remove summer months (no SSWs in this season)
def spaghetti(U_strat, start_year, end_year):
	U_strat_SOND = U_strat.extract(iris.Constraint(month = ['Sep', 'Oct', 'Nov', 'Dec']))
	U_strat_JFMAM = U_strat.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr', 'May']))
	#combine months Sep-May
	U_spaghetti = []
	for i in range(start_year, end_year):
		print i
		dummy = iris.cube.CubeList([U_strat_SOND.extract(iris.Constraint(year = i)),U_strat_JFMAM.extract(iris.Constraint(year = i+1))])
		U_spaghetti.append(dummy.concatenate_cube())
	#time array (9 months)
	time = np.linspace(0,9,len(U_spaghetti[1].data))
	print U_spaghetti
	#find standard deviation and mean for each time point
	dummy=np.zeros(len(U_spaghetti))
	std = np.zeros(len(time))
	mean = np.zeros(len(time))
	for i in range(len(U_spaghetti[1].data)):
		for j in range(1, len(U_spaghetti)-1, 1):
			dummy[j]= U_spaghetti[j].data[i]
		std[i] = np.std(dummy)
		mean[i] = np.mean(dummy)
	return mean, std, U_spaghetti



def alt_time(U_60N, start_year, end_year):
	U_SOND = U_60N.extract(iris.Constraint(month = ['Sep', 'Oct', 'Nov', 'Dec']))
	U_JFMAM = U_60N.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr', 'May']))

	#combine months Sep-May
	U_alt_time = []
	for i in range(start_year, end_year):
		dummy = iris.cube.CubeList([U_SOND.extract(iris.Constraint(year = i)),U_JFMAM.extract(iris.Constraint(year = i+1))])
		U_spaghetti.append(dummy.concatenate_cube())
	
	print U_alt_time
	

### vortex moment analysis
def vortex_calc(gph_nh, lats_nh, lons, xypoints):
	moments = vor_fast.calc_moments(gph_nh, lats_nh, lons, xypoints, hemisphere='NH',field_type='GPH', edge=3.02e4, resolution='low')	
	#moments = vor.calc_moments(gph,lats,lons,'NH','GPH',3.02e4)
	#aspect =  moments['aspect_ratio']
    	#latcent = moments['centroid_latitude']
    	kurtosis = moments['kurtosis']
    	#obj_area = moments['objective_area']

	
	return kurtosis#, aspect, latcent, obj_area


def month_warming_count(SSW_GC3, SSW_year_GC3, mon, years_GC3):
	year_count = np.empty(0)
	SSW_year_mon = SSW_year_GC3
	for i in range(len(SSW_GC3)):
		if SSW_GC3[i] != mon:
			SSW_year_mon[i] = 0
	SSW_year_mon = filter(lambda a: a != 0, SSW_year_mon)
	
	count = Counter(SSW_year_mon)
	for i in range(len(years_GC3)):
		year_count = np.append(year_count, count[str(years_GC3[i])])

	fig = plt.figure(figsize=(18, 9))
	plt.bar(years_GC3, year_count)
	plt.xlim(years_GC3[0], years_GC3[-1])
	plt.ylim(-0.1,5)
	plt.yticks([0,1,2,3,4,5])
	plt.title('Annual ' + mon + ' SSW count, GC3')
	plt.ylabel('number of SSWs')
	plt.xlabel('year')
	plt.show()
	fig.savefig('../figures/GC3_annual_' + mon + '_SSW.png')
	return year_count


def DJF_warmings_per_decade(SSW, SSW_year, years):
	SSW_year2 = SSW_year[:]
	year_count = np.empty(0)
	running = np.empty(0)

	#### count number per decade
	SSW_year_counter = Counter(SSW_year2)
	for i in range(len(years)):
		year_count = np.append(year_count, SSW_year_counter[years[i]])
	
	mean_perseason = np.mean(year_count)
	ER_perseason = np.std(year_count)/np.sqrt(len(year_count))
	#print running
	return mean_perseason, ER_perseason


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
	for i in range(len(years)-10):
		running = np.append(running, np.sum(year_count[i:i+10]))
	mean_perdec = np.mean(running)
	ER_perdec = np.std(running)
	return mean_perdec, ER_perdec
















