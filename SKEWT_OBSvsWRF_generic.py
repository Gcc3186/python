from datetime import datetime

import matplotlib.pyplot as plt
from metpy.plots import SkewT
from metpy.units import pandas_dataframe_to_unit_arrays, units
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir
from netCDF4 import Dataset
from types import SimpleNamespace

# Define a funciton that convert m/s to knot
def ms_to_kt(ms):
    return  1.944 * ms

#WRF
import wrf
from wrf import getvar, ALL_TIMES, extract_vars, ll_to_xy, to_np

dt = datetime(2013, 9, 18, 12)
# BRO
station = 'BRO' 
station_id = '72250'
lat = 25.91
lon = -97.41

## CRP
#station = 'CRP'
#station_id = '72251'#CRP
#lat = '27.76' #CRP
#lon = '-97.50'#CRP

# Read remote sounding data based on time (dt) and station
df = WyomingUpperAir.request_data(dt, station)

# Create dictionary of united arrays
data = pandas_dataframe_to_unit_arrays(df)

# Isolate united arrays from dictionary to individual variables
p = data['pressure']
T = data['temperature']
Td = data['dewpoint']
u = data['u_wind']
v = data['v_wind']
#print(v)

# Change default to be better for skew-T
fig = plt.figure(figsize=(9, 11))
#fig = plt.figure()
# Initiate the skew-T plot type from MetPy class loaded earlier
skew = SkewT(fig, rotation=45)


# Open WRF file
#wrf_filenames = ["/glade/scratch/cuchiara/wrfout_d03_2013-09-18_12:00:00"]
#wrf_filenames = ["/glade/scratch/cuchiara/wrfout_d03_2013-09-18_13:00:00"]
wrf_filenames = ["wrchem_on_12UTC"]

timex = 0 #00
#timex = 144 #12

# Open wrfouts
wrfin = [Dataset(x) for x in wrf_filenames]
# Define variables
#vars = ("tc","td","ua","va","pressure")
#for var in vars:
#    v = getvar(wrfin, var, timeidx=timex)

tc = wrf.getvar(wrfin, "tc", timeidx=timex)
td = wrf.getvar(wrfin, "td", timeidx=timex)
ua = wrf.getvar(wrfin, "ua", timeidx=timex)
va = wrf.getvar(wrfin, "va", timeidx=timex)
pressure = wrf.getvar(wrfin, "pressure", timeidx=timex)

# Find coordinate location on domain
sonde_loc = wrf.ll_to_xy(wrfin, lat, lon, timeidx=timex, squeeze=True, meta=True, stagger=None, as_int=True)

jj = sonde_loc[0]-1 
ii = sonde_loc[1]-1  

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p[:], T[:], 'black',label='Observation')
skew.plot(p[:], Td[:], 'black',linestyle='--',label='')
skew.plot_barbs(p[::5], u[::5], v[::5], y_clip_radius=0.03)


# plot WRF on Skew-T diagram
skew.plot(pressure[::3,ii,jj], td[::3,ii,jj], 'r',linestyle='--',label='')
skew.plot(pressure[::3,ii,jj], tc[::3,ii,jj], 'r',label='WRF')
skew.plot_barbs(to_np(pressure[::5,ii,jj]),to_np(ms_to_kt(ua[::5,ii,jj])),to_np(ms_to_kt(va[::5,ii,jj])),y_clip_radius=0.03,color='r')
plt.ylabel('values')


# Set some appropriate axes limits for x and y
skew.ax.set_xlim(-30, 40)
skew.ax.set_ylim(1020, 100)

# Add the relevant special lines to plot throughout the figure
skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
                       alpha=0.25, color='orangered')
skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
                         alpha=0.25, color='tab:green')
skew.plot_mixing_lines(p=np.arange(1000, 99, -20) * units.hPa,
                       linestyle='dotted', color='tab:blue')

# Add some descriptive titles
plt.title('{} Sounding'.format(station)+' ('+format(lat)+';'+format(lon)+')', loc='left')
plt.title('{}'.format(dt), loc='center')
plt.ylabel('(hPa)')
plt.xlabel('Temperature ($^\circ$C)')

# Add legend
plt.legend(loc='lower left')

# Either save or plot the figure
#plt.show()
plt.savefig('SkewT_'+str(station)+'.png')
