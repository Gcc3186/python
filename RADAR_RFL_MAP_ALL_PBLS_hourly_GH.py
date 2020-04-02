import numpy
import numpy as np
from matplotlib import pyplot as plt, patches
from matplotlib.cm import get_cmap
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
from wrf import getvar, to_np, get_cartopy, latlon_coords, cartopy_xlim,get_basemap, cartopy_ylim, Constants, interplevel, vertcross, vinterp, ALL_TIMES
from matplotlib import cm
from matplotlib.colors import from_levels_and_colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import netCDF4 as netcdf
import cartopy
import os
import glob

from netCDF4 import Dataset
from xarray import DataArray
from mpl_toolkits.basemap import Basemap

from wrf import (getvar, interplevel, vertcross,
                 vinterp,CoordPair)
from cartopy.io.shapereader import Reader

pbls = '../wrfout_links/wrchem_WRFchemV4_FINAL90lev_off_NARR_d03_2013-09-18_19:00:00'
hh = int(pbls[68:70])
print('AA '+str(hh))
ncfile1 = Dataset(str(pbls),'r')

latf = ncfile1.variables['XLAT']
lonf = ncfile1.variables['XLONG']
nlats = latf.shape[1]
nlons = latf.shape[2]

lat = latf
lon = lonf

# Create a variable "timeb" with hour:minute of output ex. t_one=12 -> timeb =['12:00','12:10'...]
t_one = hh
t_onef = hh + 2
timeb = {}
idx=0
for hour in range(t_one,t_onef):
   for minute in range(0,6):
      timeb[idx] = str(hour)+':'+str(minute)+'0'
      idx = idx + 1

ass = 0
it = 0
while it < 6:
	ass = it + 1 
	print('Working on time '+str(timeb[ass]))
	var = getvar(ncfile1,"mdbz",timeidx=it)

	# Create the figure
	fig = plt.figure(figsize=(12,9))

	ax = plt.axes()
	
	# Download and add the states and coastlines
	# See the cartopy documentation for more on this.
	states = NaturalEarthFeature(category='cultural',
                                 scale='50m',
                                 facecolor='none',
                                 name='admin_4_states_provinces_shp')

	dbz_levels = np.arange(5., 75., 5.)

	viridis = cm.get_cmap('viridis', 256)
	newcolors = viridis(np.linspace(0, 1, 256))
	dbz_rgb = np.array([[4,233,231],
                   [1,159,244], [3,0,244],
                    [2,253,2], [1,197,1],
                    [0,142,0], [253,248,2],
                    [229,188,0], [253,149,0],
                    [253,0,0], [212,0,0],
                    [188,0,0],[248,0,253],
                    [152,84,198]], np.float32)/ 255.0

	dbz_map, dbz_norm = from_levels_and_colors(dbz_levels, dbz_rgb,
                                           extend="max")

	newcolors = dbz_rgb
	newcmp = ListedColormap(newcolors)

	levels = np.arange(5, 75, 5.)

	# Convert the lat/lon coordinates to x/y coordinates in theprojectionspace
	x, y = to_np(lonf), to_np(latf)

	if it < 101:
                #Define map boundaries
		map = Basemap(llcrnrlon=-99.0,llcrnrlat=26.5,urcrnrlon=-95.0,urcrnrlat=29.5,	resolution = 'h')
	

	plt.contourf(lonf[0,:,:], latf[0,:,:], var[:,:], var[:,:],levels=levels,cmap=dbz_map)

	map.drawcoastlines()
	map.drawparallels(np.arange(20.5,39,.5),labels=[1,0,0,0],linewidth=0,fontsize=18)
	map.drawmeridians(np.arange(0,360,1),labels=[0,0,0,1],linewidth=0,fontsize=18)
	cbar = map.colorbar(location='right')
	cbar.set_label('(dBZ)',rotation=270,fontsize=12)
	cbar.ax.tick_params(labelsize=12) 
	plt.title('Maximun Reflectivity - '+str(timeb[ass]+' (UTC)'),fontsize=18)
	# Shape file with country boundaries (need to download)
	#Plot states
	#map.readshapefile('/glade/work/cuchiara/ncl_wks/Sep02/shp/USA_adm2', 'comarques')
	# Plot roads
	#map.readshapefile('/glade/work/cuchiara/ncl_wks/Sep02/shp/major_us_roads/tl_2016_us_primaryroads', 'comarques2',color='grey',linewidth=.5)
	#plt.clabel(levels, inline=1, fontsize=10, fmt="%i")
	#plt.show()
	plt.savefig('RFL_'+str(timeb[ass])+'.png')
	plt.clf()
	plt.close()
	del(var)
	ass = ass + 1
	it += 1
