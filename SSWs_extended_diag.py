import iris
import iris.coord_categorisation as coord_cat
import matplotlib.pyplot as plt
from matplotlib import gridspec
from GC3_diagnostics_functions import *
import numpy as np


def callback(cube, field, filename): 
    # Remove the history attribute. 
    del cube.attributes['history']
    del cube.attributes['valid_max']
    del cube.attributes['valid_min']

def SSW_counter(U_strat, start_year, end_year, thresh):

	#identify SSWs from U_strat cube. defition of an SSW: U becomes negative at 60N 10hPa in months NDJFM (20 days after an SSW, another can occur)
	#restrict months
	SSW = iris.cube.CubeList([])
	U_SSW_seasons, SSW_month, final, indexi, indexj, SSW_year = [], [], [], [], [], []
	U_strat_ND = U_strat.extract(iris.Constraint(month = ['Nov', 'Dec']))
	U_strat_JFMA = U_strat.extract(iris.Constraint(month = ['Jan', 'Feb', 'Mar', 'Apr']))	
	U_SSW_seasons.append(U_strat_JFMA.extract(iris.Constraint(year = start_year)))

	for i in range(start_year, end_year):
		print i
		dummy2 = iris.cube.CubeList([U_strat_ND.extract(iris.Constraint(year = i)),U_strat_JFMA.extract(iris.Constraint(year = i+1))])
		U_SSW_seasons.append(dummy2.concatenate_cube())
	
	for j in range(len(U_SSW_seasons)):
		print j	
		i=0
		while i < len(U_SSW_seasons[j].data):
			#print i,j
			if U_SSW_seasons[j].data[i] < thresh and U_SSW_seasons[j].data[i-1] >= thresh:
				print "warming detected"
				SSW_month.append(str(U_SSW_seasons[j].coord('month').points[i]))
				SSW_year.append(U_SSW_seasons[j].coord('season_year').points[i])
				SSW.append(U_SSW_seasons[j][i])
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
		
		if final != []:	
			for k in range(final[-1],len(U_SSW_seasons[j].data)-9):
				if np.all(U_SSW_seasons[j].data[k:k+9] > thresh):
					break
				elif k== len(U_SSW_seasons[j].data)-10: 
					#print 'final warming found'
					del SSW[-1]
					del SSW_month[-1]
					del SSW_year[-1]
					break
				else:
					k=k	
		else:	
			j = j
			
			
	for i in range(len(SSW)):
		if SSW_month[i] == 'Apr':
			SSW_year[i] = 0
			SSW[i] = 0
	SSW_year = filter(lambda a: a != 0, SSW_year)
    	SSW_month = filter(lambda a: a != 'Apr', SSW_month)
	SSW = filter(lambda a: a != 0, SSW)

	print SSW_year, len(SSW_year)
	return SSW, SSW_year, SSW_month


def filter_nov_years(SSW, SSW_month, SSW_year):

	for i in range(len(SSW)):
		print i
		if SSW_month[i] == 'Nov':
			SSW_month[i] = 0
			SSW[i] = 0
			if SSW_year[i+1] == SSW_year[i]:
				SSW_year[i+1] = 0
				SSW_year[i] = 0
				SSW_month[i+1] = 0
				SSW[i+1] = 0
			
			elif SSW_year[i+1] == SSW_year[i] + 1 and SSW_month[i+1] == 'Jan' or 'Feb' or 'Mar':
				SSW_year[i+1] = 0
				SSW_year[i] = 0
				SSW_month[i+1] = 0
				SSW[i+1] = 0

		if SSW_month[i] == 'Mar':
			SSW_month[i] = 0
			SSW[i] = 0
			SSW_year[i] = 0
				
	SSW_year = filter(lambda a: a != 0, SSW_year)
    	SSW_month = filter(lambda a: a != 0, SSW_month)
	SSW = filter(lambda a: a != 0, SSW)

	return SSW, SSW_year, SSW_month			



def DJF_warmings_per_decade(SSW, SSW_year, years):
	SSW_year2 = SSW_year[:]
	year_count = np.empty(0)
	running = np.empty(0)

	#### count number per decade
	SSW_year_counter = Counter(SSW_year2)
	print SSW_year_counter
	for i in range(len(years)):
		year_count = np.append(year_count, SSW_year_counter[years[i]])
	print len(years)
	for i in range(len(years)-10):
		running = np.append(running, np.sum(year_count[i:i+10]))
	mean_perdec = np.mean(running)
	ER_perdec = np.std(running)
	#print running
	return year_count



def warmings_per_winter(year_count):
	per_season = np.mean(year_count)
	Er = np.std(year_count)/np.sqrt(len(year_count))
	
	return per_season, Er	
	

def monte_carlo(year_count, years):
	SSW_perseason_rean = np.empty(0)
	for i in range(len(years)-rean_per):
		print i
		start = i
		period = year_count[start:start+rean_per]
		SSW_perseason_rean = np.append(SSW_perseason_rean, np.mean(period))
	return SSW_perseason_rean




fname = '/gws/nopw/j04/aopp/odimdore/long_runs/CMIP_data/U/N96/daily/u-ar766_U_zonmeans_1850-2350.nc'
cubes = iris.load(fname, lat_constraint)#, callback=callback)
U_N96 = add_times(cubes[0], 'time')
U_N96 = U_N96[:,15,0]
start_yearN96 = U_N96.coord('year').points[0]
end_yearN96 = U_N96.coord('year').points[-1]
years_N96 = np.arange(start_yearN96, end_yearN96+1, 1)

fname = '/gws/nopw/j04/aopp/odimdore/long_runs/CMIP_data/U/N216/daily/total.nc'
cubes = iris.load(fname, lat_constraint)# callback=callback)
U_N216 = add_times(cubes[0], 'time')
U_N216 = U_N216[:,0,0]
start_yearN216 = U_N216.coord('year').points[0]
end_yearN216 = U_N216.coord('year').points[-1]
years_N216 = np.arange(start_yearN216, end_yearN216+1, 1)

fname = '/gws/nopw/j04/aopp/strat_clim_data/CMIP6/UKESM1/picontrol/u_pl/daily/UKESM*'
cubes = iris.load(fname, lat_constraint, callback=callback)
del cubes[3].attributes['nco_openmp_thread_number']
del cubes[3].attributes['NCO']
U_UKESM = cubes.concatenate_cube()
U_UKESM = add_times(U_UKESM, 't')
U_UKESM = U_UKESM[:,0,0]
start_yearUKESM = U_UKESM.coord('year').points[0]
end_yearUKESM = U_UKESM.coord('year').points[-1]
years_UKESM = np.arange(start_yearUKESM, end_yearUKESM+1, 1)


SSW_ERA = ['Feb', 'Feb', 'Mar', 'Dec', 'Feb', 'Jan', 'Jan', 'Dec', 'Mar', 'Feb', 'Dec', 'Feb', 'Mar', 'Feb', 'Dec', 'Jan', 'Jan', 'Jan', 'Feb', 'Feb', 'Jan', 'Feb', 'Mar', 'Jan']
SSW_year_ERA = [1979, 1980, 1981, 1981, 1984, 1985, 1987, 1987, 1988, 1989, 1998, 1999, 2000, 2001, 2001, 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2010, 2013]
SSW_ERA40 = ['Jan', 'Dec', 'Feb', 'Jan', 'Nov', 'Mar', 'Jan', 'Jan', 'Mar', 'Jan', 'Jan']
SSW_year_ERA40 = [1963, 1965, 1966, 1968, 1968, 1969, 1970, 1971, 1971, 1973, 1977]
SSW_ERA2 = ['Feb', 'Feb', 'Mar']
SSW_year_ERA2 = [2017, 2018, 2018]

###combine reanalysis datests
SSW_ERA = SSW_ERA40 + SSW_ERA + SSW_ERA2
SSW_year_ERA = SSW_year_ERA40 + SSW_year_ERA + SSW_year_ERA2

start_yearERA = 1959 
end_yearERA = 2018
years_ERA = np.arange(start_yearERA, end_yearERA, 1)
rean_per = len(years_ERA)


#ERA
test1, djfSSW_year_ERA, test2= filter_nov_years(SSW_ERA, SSW_ERA, SSW_year_ERA)
djfSSW_ERA = filter(lambda a: a != 'Mar', SSW_ERA)
djfSSW_ERA = filter(lambda a: a != 'Nov', djfSSW_ERA)


SSW_ERA_djf, ERA_year_djf, ERA_month_djf= filter_nov_years(SSW_ERA, SSW_ERA, SSW_year_ERA)
year_count_ERA = DJF_warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA)
ERA_per_season, ErERA = warmings_per_winter(year_count_ERA)

djfyear_count_ERA = DJF_warmings_per_decade(SSW_ERA_djf, ERA_year_djf, years_ERA)
djfERA_per_season, djfErERA = warmings_per_winter(djfyear_count_ERA)




#N96
SSW_N96, SSW_year_N96, SSW_month_N96 = SSW_counter(U_N96, start_yearN96, end_yearN96, 0)
#djf months
SSW_N96_djf, N96_year_djf, N96_month_djf= filter_nov_years(SSW_N96, SSW_month_N96, SSW_year_N96)
djfyear_count_N96 = DJF_warmings_per_decade(N96_month_djf, N96_year_djf, years_N96)
djfper_wint_N96, djfErN96 = warmings_per_winter(djfyear_count_N96)
djfSSW_perseason_reanN96 = monte_carlo(djfyear_count_N96, years_N96)
#all months
year_count_N96 = DJF_warmings_per_decade(SSW_month_N96, SSW_year_N96, years_N96)
per_wint_N96, ErN96 = warmings_per_winter(year_count_N96)
SSW_perseason_reanN96 = monte_carlo(year_count_N96, years_N96)


#N216
SSW_N216, SSW_year_N216, SSW_month_N216 = SSW_counter(U_N216, start_yearN216, end_yearN216, 0)
#djf months
SSW_N216_djf, N216_year_djf, N216_month_djf= filter_nov_years(SSW_N216, SSW_month_N216, SSW_year_N216)
djfyear_count_N216 = DJF_warmings_per_decade(N216_month_djf, N216_year_djf, years_N216)
djfper_wint_N216, djfErN216 = warmings_per_winter(djfyear_count_N216)
djfSSW_perseason_reanN216 = monte_carlo(djfyear_count_N216, years_N216)
#all months
year_count_N216 = DJF_warmings_per_decade(SSW_month_N216, SSW_year_N216, years_N216)
per_wint_N216, ErN216 = warmings_per_winter(year_count_N216)
SSW_perseason_reanN216 = monte_carlo(year_count_N216, years_N216)


#UKESM
SSW_UKESM, SSW_year_UKESM, SSW_month_UKESM = SSW_counter(U_UKESM, start_yearUKESM, end_yearUKESM, 0)
#djf months
SSW_UKESM_djf, UKESM_year_djf, UKESM_month_djf= filter_nov_years(SSW_UKESM, SSW_month_UKESM, SSW_year_UKESM)
djfyear_count_UKESM = DJF_warmings_per_decade(UKESM_month_djf, UKESM_year_djf, years_UKESM)
djfper_wint_UKESM, djfErUKESM = warmings_per_winter(djfyear_count_UKESM)
djfSSW_perseason_reanUKESM = monte_carlo(djfyear_count_UKESM, years_UKESM)
#all months
year_count_UKESM = DJF_warmings_per_decade(SSW_month_UKESM, SSW_year_UKESM, years_UKESM)
per_wint_UKESM, ErUKESM = warmings_per_winter(year_count_UKESM)
SSW_perseason_reanUKESM = monte_carlo(year_count_UKESM, years_UKESM)





series = np.array([per_wint_N96, per_wint_UKESM, ERA_per_season])
Er_series = np.array([ErN96, ErUKESM, ErERA])
index = np.linspace(0, 2, 3)

fig, ax = plt.subplots(figsize = (5,5), dpi=300)
bar_width = 0.35
plt.bar(index, series, bar_width, yerr=Er_series,  capsize = 5, ecolor = 'black', color = ('r','b','cyan'), linewidth = 2.0)
ax.set_xticks(index)
ax.set_xticklabels((['GC3 N96','UKESM','ERA']), fontsize = 10)
plt.title('Mean NDJFM SSWs per winter', fontsize = 15)
plt.ylabel('SSWs per winter', fontsize = 15)
plt.yticks(np.arange(0,1,0.1))
plt.ylim(0,1)
plt.show()
fig.savefig('../../figures/SSWs_per_winter_bar_for_ACSIS.png')





binwidth = 0.017
fig, ax = plt.subplots(dpi=300)
gs = gridspec.GridSpec(3, 2)
#fig.title('DJF warmings per winter', fontsize = 25)
plt1 = plt.subplot(gs[0,0])
plt.hist(SSW_perseason_reanN96, bins=np.arange(min(SSW_perseason_reanN96), max(SSW_perseason_reanN96) + binwidth, binwidth), color = 'r', alpha = 0.4, density = True, label = 'GC3 N96')
bins=np.arange(min(SSW_perseason_reanN96), max(SSW_perseason_reanN96) + binwidth, binwidth)
#(mu, sigma) = norm.fit(SSW_perseason_reanN96)
#y = mlab.normpdf( bins, mu, sigma)
#plt.plot(bins, y, color = 'r')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black', linestyle = '--', linewidth = 0.5)
plt.ylabel('frequency (%)', fontsize = 5)
plt.xlim(0.25,1)
plt.yticks(np.arange(8), fontsize = 5)
plt.xticks(fontsize = 5)
plt.ylim(0,8)
plt.legend(fontsize = 5)

plt3 = plt.subplot(gs[1,0])
plt.hist(SSW_perseason_reanN216, bins=np.arange(min(SSW_perseason_reanN216), max(SSW_perseason_reanN216) + binwidth, binwidth), color = 'g', alpha = 0.4, density = True, label = 'GC3 N216')
bins=np.arange(min(SSW_perseason_reanN216), max(SSW_perseason_reanN216) + binwidth, binwidth)
#(mu, sigma) = norm.fit(SSW_perseason_reanN216)
#y = mlab.normpdf( bins, mu, sigma)
#plt.plot(bins, y, color = 'g')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black', linestyle = '--', linewidth = 0.5)
plt.ylabel('frequency (%)', fontsize = 5)
plt.xlim(0.25,1)
plt.yticks(np.arange(8), fontsize = 5)
plt.xticks(fontsize = 5)
plt.ylim(0,8)
plt.legend(fontsize = 5)


plt5 = plt.subplot(gs[2,0])
plt.hist(SSW_perseason_reanUKESM, bins=np.arange(min(SSW_perseason_reanUKESM), max(SSW_perseason_reanUKESM) + binwidth, binwidth), color = 'b', alpha = 0.4, density = True, label = 'UKESM1 N96')
bins=np.arange(min(SSW_perseason_reanUKESM), max(SSW_perseason_reanUKESM) + binwidth, binwidth)
#(mu, sigma) = norm.fit(SSW_perseason_reanUKESM)
#y = mlab.normpdf(bins, mu, sigma)
#plt.plot(bins, y, color = 'b')
plt.plot([ERA_per_season, ERA_per_season],[0,10], label = 'SSWs in ERA', color = 'black', linestyle = '--', linewidth = 0.5)
plt.xlabel('SSWs per season', fontsize = 5)
plt.ylabel('frequency (%)', fontsize = 5)
plt.xlim(0.25,1)
plt.yticks(np.arange(8), fontsize = 5)
plt.xticks(fontsize = 5)
plt.ylim(0,8)
plt.legend(fontsize = 5)
#plt.tight_layout()
plt.show()
fig.savefig('../../figures/SSWs_perseason_RA_period_for_transfer.png')








months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']


def mon_warmings_per_decade(SSW, SSW_year, years, mon):
	SSW_year2 = SSW_year[:]
	year_count = np.empty(0)
	running = np.empty(0)
	#### filter to warmings of a single month
	for i in range(len(SSW)):
		if SSW[i] != str(mon):
			SSW_year2[i] = 0
	SSW_year2 = filter(lambda a: a != 0, SSW_year2)
	
	#### count number per season
	SSW_year_counter = Counter(SSW_year2)
	for i in range(len(years)):
		year_count = np.append(year_count, SSW_year_counter[years[i]])
	
	per_season = np.mean(year_count)
	Er = np.std(year_count)/np.sqrt(len(year_count))
	return per_season, Er
	
	


N96, N216, ERA, UKESM, ErN96, ErN216, ErERA, ErUKESM = [], [], [], [], [], [], [], []
for i in range(len(months)):
	N96.append(mon_warmings_per_decade(SSW_N96, SSW_year_N96, years_N96, months[i])[0])
	UKESM.append(mon_warmings_per_decade(SSW_UKESM, SSW_year_UKESM, years_UKESM, months[i])[0])
	#N216.append(mon_warmings_per_decade(SSW_N216, SSW_year_N216, years_N216, months[i])[0])
	ERA.append(mon_warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA, months[i])[0])
	ErN96.append(mon_warmings_per_decade(SSW_month_N96, SSW_year_N96, years_N96, months[i])[1])
	#ErN216.append(mon_warmings_per_decade(SSW_month_N216, SSW_year_N216, years_N216, months[i])[1])
	ErERA.append(mon_warmings_per_decade(SSW_ERA, SSW_year_ERA, years_ERA, months[i])[1])
	ErUKESM.append(mon_warmings_per_decade(SSW_month_UKESM, SSW_year_UKESM, years_UKESM, months[i])[1])


#SSW_N216_counter = Counter(SSW_month_N216)
SSW_N96_counter = Counter(SSW_month_N96)
SSW_UKESM_counter = Counter(SSW_month_UKESM)
SSW_ERA_counter = Counter(SSW_ERA)

SSWs_ERA = np.array(([SSW_ERA_counter['Nov'],SSW_ERA_counter['Dec'],SSW_ERA_counter['Jan'],SSW_ERA_counter['Feb'],SSW_ERA_counter['Mar']])).astype(np.float)/len(years_ERA)
SSWs_N96 = np.array(([SSW_N96_counter['Nov'],SSW_N96_counter['Dec'],SSW_N96_counter['Jan'],SSW_N96_counter['Feb'],SSW_N96_counter['Mar']])).astype(np.float)/len(years_N96)
#SSWs_N216 = np.array(([SSW_N216_counter['Nov'],SSW_N216_counter['Dec'],SSW_N216_counter['Jan'],SSW_N216_counter['Feb'],SSW_N216_counter['Mar']])).astype(np.float)/len(years_N216)
SSWs_UKESM = np.array(([SSW_UKESM_counter['Nov'],SSW_UKESM_counter['Dec'],SSW_UKESM_counter['Jan'],SSW_UKESM_counter['Feb'],SSW_UKESM_counter['Mar']])).astype(np.float)/len(years_UKESM)

index = np.linspace(0, 10, 5)
fig, ax = plt.subplots(dpi=300)
bar_width = 0.35

warm1 = ax.bar(index - bar_width, SSWs_N96, bar_width, yerr=ErN96,  capsize = 5, ecolor = 'black', color = 'r', linewidth = 2.0)
#warm2 = ax.bar(index - 0.5*bar_width, SSWs_N216, bar_width, yerr=ErN216, capsize = 5, ecolor = 'black', color = 'g', linewidth = 2.0)
warm3 = ax.bar(index, SSWs_ERA, bar_width, yerr=ErERA, capsize = 5, ecolor = 'black', color = 'cyan', linewidth = 2.0)
warm4 = ax.bar(index + bar_width, SSWs_UKESM, bar_width, yerr=ErUKESM, capsize = 5, ecolor = 'black', color = 'b', linewidth = 2.0)

ax.set_xticks(index)
ax.set_xticklabels(('Nov', 'Dec', 'Jan', 'Feb', 'Mar'))
ax.legend((warm1[0],warm3[0],warm4[0]),('GC3 N96', 'ERA', 'UKESM'))

plt.title('Mean SSWs per winter', fontsize = 15)
plt.ylabel('SSWs per winter', fontsize = 15) 
fig.tight_layout()
plt.show()
fig.savefig('../../figures/month_hist_for_ACSIS.png')

