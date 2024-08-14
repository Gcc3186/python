# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %%
from __future__ import (absolute_import, division, print_function)
from mpl_toolkits.basemap import Basemap, shiftgrid # Requires python 3.7.6
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
plt.rcParams['text.usetex'] = False
from PseudoNetCDF.camxfiles.Memmaps import uamiv
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas
# %%
# %%
f = netCDF4.Dataset('/npsnas7/NEIC_2016_V1/simulations_output/CAMx.12US2.CAMx.v7.00.MPICH3.NCF4.pgfomp.20151222.avrg.grd01.nc')

netCDF4.Dataset('/npsnas7/NEIC_2016_V1/simulations_output/CAMx.12US2.CAMx.v7.00.MPICH3.NCF4.pgfomp.20151222.avrg.grd01.nc')

# Get Ozone Variable
latf = f.variables['latitude']
lonf = f.variables['longitude']
nlats = latf.shape[0]
nlons = lonf.shape[0]
#lat = data.createDimension("lat", nlats)###
lat = latf
#lon = data.createDimension("lon", nlons)
lon = lonf
# %%
import glob
for filepath in glob.iglob('/npsnas7/NEIC_2016_V1/simulations_output/CAMx.12US2.CAMx.v7.00.MPICH3.NCF4.pgfomp.201601*.avrg.grd01'):
    #print(filepath)
    it = 0
    k = 0
    while it < 24:
        #camx = uamiv('/npsnas7/NEIC_2016_V1/simulations_output/CAMx.12US2.CAMx.v7.00.MPICH3.NCF4.pgfomp.20160101.avrg.grd01')
        camx = uamiv(filepath)
        o3 = camx.variables['O3'][it,0,:,:]*1000
        vvm = 75
        
        m = Basemap(
        projection='lcc',resolution = 'h',
        lon_0=f.P_GAM, lat_0=f.YCENT, lat_1=f.P_ALP, lat_2=f.P_BET,
        llcrnrx=f.XORIG, llcrnry=f.YORIG,
        urcrnrx=f.XORIG + f.NCOLS * f.XCELL,
        urcrnry=f.YORIG + f.NROWS * f.YCELL)

        x = np.arange(f.NCOLS + 1) * f.XCELL
        y = np.arange(f.NROWS + 1) * f.YCELL


        # create the figure.
        fig=plt.figure(figsize=(15,20))

        # add an axes.
        ax = fig.add_axes([0.1,0.1,0.9,0.9])
        ax.set_facecolor('lightgrey')
        # associate this axes with the Basemap instance.
        m.ax = ax

        # plot tile plot with pcolor
        
        from matplotlib.colors import ListedColormap
        WhGrYlBu = ListedColormap(['#ffffff', '#b7f6ff', '#70edff', '#29e4ff', '#00e1fb', '#0fffc6', '#3bffa4', '#68ff82', '#94ff60', '#c0ff3e', '#edff1c', '#fff400', '#ffc700', '#ff9b00', '#ff6e00', '#ff4200', '#ff1500', '#e80000', '#bb0000', '#8f0000'])
        #levels = np.arange(0, 135, 5.)
        #toplot = np.ma.masked_values(o3[:,:]*1000)
        plt.pcolormesh(x, y, o3[:,:],cmap = WhGrYlBu,vmin=0,vmax=vvm)
        #plt.contourf(lonf[:,:], latf[:,:], o3[:,:]*1000,levels=levels,cmap = WhGrYlBu,extend="both")
        #plt.contourf(x, y, o3[:,:]*1000,levels=levels,cmap = WhGrYlBu,extend="both")
    
    
        # Add a colorbar using proportional spacing, but
        #cbar.ax.xaxis.get_offset_text().set_position((1,100))
        cb = m.colorbar(location='right')
        cb.set_label("(ppbv)",fontsize=18)
        cb.ax.tick_params(labelsize=18)
        #cb.set_ticks([0,5])
    
        # Add units to the colorbar
        # Sites file info
        filename = "/npsnas7/NEIC_2016_V1/Observations/SITEO3.csv"
        colnames = ['Site', 'Code', 'Dataset', 'Country', 'State', 'County', 'AQSCode', 'Latitude', 'Longitude', 'Elevation', 'DemographicCode', 'LandUseCode']
        data = pandas.read_csv(filename, names=colnames)
        s = data['Site']
        s = s.array
        c = data['Code']
        c = c.array
        lon = data['Longitude']
        lon = lon.array
        lat = data['Latitude']
        lat = lat.to_numpy(dtype='object')
        aqs = data['AQSCode']
        aqs = aqs.to_numpy(dtype='object')
            
        # Observation file
        filevar = "/npsnas7/NEIC_2016_V1/Observations/O3.csv"
        colname = ['Site', 'Code', 'Date', 'Var']
        data_var = pandas.read_csv(filevar, names=colname)

        sc = data_var['Site']
        sc = sc.to_numpy(dtype='object')
        d = data_var['Date']
        d = d.array
        var = data_var['Var']
        var = var.array
        
        cx = d.shape

        lara = [0] * cx[0]
        lora = [0] * cx[0]
        vara = [0] * cx[0]
        a = 0
        
        for i, j in enumerate(d):
            if d[i] == '1/1/2016 '+str(k)+':00':
                #print(str('1/1/2016 '+str(i)+':00'))
                #print(sc[i])
                ind = int(np.argwhere(aqs == float(sc[i])))
                lara[a] = float(lat[ind])
                lora[a] = float(lon[ind])
                vara[a] = float(var[i])
                a = a + 1
                #print(d[i]+' '+str(var[i])+' '+str(sc[i])+' '+str(c[ind])+' '+str(lat[ind])+' '+str(lon[ind]))
        k = k + 1
        
        # plot blue dot on Houston, Baton Rouge, and Atlanta
        
        #def add_dot(lono, lato, valo):
        #    xpt,ypt = m(lono,lato) 
        #    #m.plot([xpt],[ypt],'o')
        #    m.scatter([xpt],[ypt],s=50,c=[valo],cmap=WhGrYlBu,edgecolors='black',vmin=0,vmax=vvm)
        #    ax.annotate(label, xy=(xpt, ypt), xytext=(xpt+1e5, ypt+1e5),
        #    #bbox=dict(boxstyle="round4", fc="w"),
        #   # arrowprops=dict(facecolor='black'),
        #    #)
        #for x in range(1000):
        #    add_dot(lora[int(x)],lara[int(x)],vara[int(x)])
        #import math
        #for x in range(799):
            #if vara[int(x)] == -999:
            #    print('Skipped!')
            #elif vara[int(x)] > 150:
            #    print('Odd!')
            #else:
                #add_dot(lora[int(x)],lara[int(x)],vara[int(x)])
                #print(vara[int(x)])
        #plt.scatter(-90.0,33.0,s=40,c=0,cmap='jet',edgecolors='black',vmin=0,vmax=200)
        
        # draw coastlines and political boundaries.
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        # draw parallels and meridians.
        # label on left, right and bottom of map.
        parallels = np.arange(20.,60,10.)
        m.drawparallels(parallels,labels=[1,0,0,1],fontsize=18, linewidth=0)
        meridians = np.arange(-120., 70.,10.)
        m.drawmeridians(meridians,labels=[1,0,0,1],fontsize=18, linewidth=0)

        #m.readshapefile('/home/cuchiara/shp/USA_adm1', 'comarques')
        # Plot roads
        #map.readshapefile('/glade/work/cuchiara/ncl_wks/Sep02/shp/major_us_roads/tl_2016_us_primaryroads', 'comarques2',color='grey',linewidth=.5)

        # handle time
        import datetime
        jday = camx.variables['TFLAG'][:,0,0]
        hour = camx.variables['TFLAG'][it,0,1]
        date = datetime.datetime.strptime(str(jday[1]), '%Y%j').strftime("%d-%m-%Y")

        # set title.
        if it < 10:
            ax.set_title('O$_3$ as predicted by the CAMx v7 - '+str(date)+' - 0'+str(int(hour))[0:1]+' (UTC)',fontsize=18)
        else:    
            ax.set_title('O$_3$ as predicted by the CAMx v7 - '+str(date)+' - '+str(int(hour))[0:2]+' (UTC)',fontsize=18)
        import textwrap
        #histstr = 'Processing: %s' % '\n'.join(textwrap.wrap(camx.history.strip(), 140))
        print(str(it))
        #fig.text(0.01, 0.01, histstr, horizontalalignment = 'left', verticalalignment = 'bottom', size = 8)
        #plt.savefig('O3_'+str(it)+'.png',bbox_inches='tight')
        #plt.close()
        plt.show()
        del(o3)
        it += 1


# %%
d


# %%
'1/1/2016 '+str(i)+':00'


# %%
from datetime import datetime

timestamp = '1/1/2016 0:00'#1230906310
date_time = datetime.fromtimestamp(timestamp)

print("Date time object:", date_time)

da = date_time.strftime("%-m/%-d/%Y, %-H:%M")
print("Output 2:", da)	


# %%



