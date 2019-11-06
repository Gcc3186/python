#!/glade/u/apps/ch/opt/python/2.7.14/intel/17.0.1/pkg-library/20180329/bin/python


import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy.stats import norm
import netCDF4 as netcdf
import os
import glob
import matplotlib.pyplot as plt


radarfiles1 = glob.glob('radar_nexrard/*17*.nc')
radarfiles2 = glob.glob('radar_nexrard/*18*.nc')
radarfiles3 = glob.glob('radar_nexrard/*19*.nc')
radarfiles4 = glob.glob('radar_nexrard/*T200*.nc')
radarfiles5 = glob.glob('radar_nexrard/*T201*.nc')
radarfiles6 = glob.glob('radar_nexrard/*T202*.nc')
radarfiles7 = glob.glob('radar_nexrard/*T203*.nc')

radarfiles = radarfiles1 + radarfiles2 + radarfiles3 + \
    radarfiles4 + radarfiles5 + radarfiles6 + radarfiles7
# radarfiles = radarfiles3

timesteps = len(radarfiles)

# Open one file to get x-y dimensions of the data
f = netcdf.Dataset(radarfiles[0], 'r')
refl = f.variables['ZH']
refl = np.asarray(refl[0, :, 150:440, 135:300])

zdim = refl.shape[0]
jdim = refl.shape[1]
idim = refl.shape[2]
# print(zdim)
# print(idim)
# print(jdim)
# quit()
f.close()


moderate_columns = list()
heavy_columns = list()
hgt_columns_heavy = list()
hgt_columns_moderate = list()
# heavy_count = 0
# moderate_count = 0

for rfile in radarfiles:

 #   print 'Opening File : ',rfile
    f = netcdf.Dataset(rfile, 'r')
    refl = f.variables['ZH']
    refl = np.asarray(refl)
    height = f.variables['Altitude']
    height = np.asarray(height)
    sl3d = f.variables['SL3D']
    sl3d = np.asarray(sl3d)
    heavy_count = 0
    moderate_count = 0

    # Now loop through grid points looking for "Moderate" and "Heavy" columns

    for jj in range(jdim):
        for ii in range(idim):

            avg_refl_1to4km = np.mean(refl[0, 0:4, jj, ii])
            if (avg_refl_1to4km >= 40.0 and (sl3d[0, jj, ii] == 1 or sl3d[0, jj, ii] == 2)):
                heavy_columns.append(refl[0, :, jj, ii])
                hgt_columns_heavy.append(height[:])
                heavy_count = heavy_count + 1
            if (avg_refl_1to4km >= 20.0 and avg_refl_1to4km <= 30.0):
                moderate_columns.append(refl[0, :, jj, ii])
                hgt_columns_moderate.append(height[:])
                moderate_count = moderate_count + 1

   # print '*******************************************'
   # print 'For time ',rfile
   # print 'Heavy columns = ',heavy_count
   # print 'Moderate columns = ',moderate_count
   # print 'Percentage of domain covered by heavy columns = ',(heavy_count/(jdim*idim))*100.0,' %'


# Conver to numpy arrays for stats
heavy_columns = np.asarray(heavy_columns)
moderate_columns = np.asarray(moderate_columns)
hgt_columns_heavy = np.asarray(hgt_columns_heavy)
hgt_columns_moderate = np.asarray(hgt_columns_moderate)


hgtrange = np.arange(0.5, 17.5, 1.0)
hgtrange_contour = np.arange(1.0, 17.0, 1.0)
dbzrange = np.arange(0, 72, 2)
dbzrange_contour = np.arange(1, 71, 2)


# Moderate CFAD
H_moderate, xedges2, yedges2 = np.histogram2d(hgt_columns_moderate.flatten(
), moderate_columns.flatten(), bins=(hgtrange, dbzrange), normed=True)

# Heavy CFAD
H_heavy, xedges3, yedges3 = np.histogram2d(hgt_columns_heavy.flatten(
), heavy_columns.flatten(), bins=(hgtrange, dbzrange), normed=True)

# print 'H_heavy shape is ',H_heavy.shape


# Now to normalize by maximum frequency at each altitude
# print 'Allocating arrays'
Nc_moderate = np.zeros((H_moderate.shape[0], H_moderate.shape[1]))
Nc_heavy = np.zeros((H_heavy.shape[0], H_heavy.shape[1]))
# print 'Done with Allocation...time to assign values'
for h in range(H_heavy.shape[0]):
    Nc_moderate[h, :] = np.amax(H_moderate[h, :])
    Nc_heavy[h, :] = np.amax(H_heavy[h, :])

print('Done assigning values to arrays')
# Now calculate percentage of maximum frequency
H_moderate = (H_moderate / Nc_moderate) * 100.
H_heavy = (H_heavy / Nc_heavy) * 100.

font = {'family': 'Arial',
        'weight': 'bold',
        'size': 22}

plt.rc('font', **font)

# Plot Moderate CFAD
f = plt.figure(figsize=(11.0, 8.5))
ax = f.add_subplot(111)
levels = np.arange(10, 110, 10)
levels50 = [50]
nlevels = np.arange(0.0, (0.0120 + (5.0 * 0.0015)), 0.0015)
CS = plt.contourf(dbzrange_contour, hgtrange_contour,
                  H_moderate, levels, cmap=plt.cm.jet, extend='min')
# CS = plt.contourf(dbzrange_contour,hgtrange_contour,H_moderate,nlevels,cmap=plt.cm.jet,extend = 'min')
CS.cmap.set_under('white')
plt.colorbar()
plt.grid(True)
plt.xlabel('Reflectivity (dBZ)', **font)
plt.ylabel('Altitude (km)', **font)
plt.savefig('Moderate_CFAD_OBS_Normalized_by_Numpy_mycolors_Cheyenne_ORIG.png')
plt.clf()


f = plt.figure(figsize=(11.0, 8.5))
ax = f.add_subplot(111)
levels = np.arange(10, 110, 10)
levels50 = [50]
CS = plt.contourf(dbzrange_contour, hgtrange_contour, H_heavy,
                  levels, cmap=plt.cm.jet, extend='min')
# CS = plt.contourf(dbzrange_contour,hgtrange_contour,H_heavy,nlevels,cmap=plt.cm.jet,extend = 'min')
CS.cmap.set_under('white')
plt.colorbar()
plt.grid(True)
plt.xlabel('Reflectivity (dBZ)', **font)
plt.ylabel('Altitude (km)', **font)
plt.savefig('Heavy_CFAD_OBS_Normalized_by_Numpy_mycolors_Cheyenne_ORIG.png')
quit()


'''

        byte SL3D(time, Latitude, Longitude) ;
                string SL3D:long_name = "SL3D Echo Classification" ;
                string SL3D:units = " " ;
                SL3D:missing_value = 0 ;
                string SL3D:values = "1 - convective updraft, 2 - convection, 3 - precipitating stratiform, 4 - non-precipitating stratiform, 5 - anvil" ;

'''
