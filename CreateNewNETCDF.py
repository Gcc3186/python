import os
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, to_np, interplevel

# Set WRF file paths
hrr = "17"
output_filename = "WRFOUT_COREA.nc"
output_dir = ""
input_dir = "wrfouts/"
filename_off = "wrfchem_rst_noscav_d03_2013-09-02_18:15:00"
filename = "wrfchem_rst_rfstd_d03_2013-09-02_18:15:00"

# Remove existing output file if present
if os.path.exists(output_dir + output_filename):
    os.remove(output_dir + output_filename)

# Open WRF files
file_off = Dataset(input_dir + filename_off, "r")
file = Dataset(input_dir + filename, "r")

# Define slicing indices
str_i, end_i = 165, 182
str_j, end_j = 337, 357
str_t, end_t = 3, 10

# Read in data
print("Importing variables from file: {}. This may take a while!".format(filename_off))

times = file.variables["Times"][str_t:end_t, :]
lat = file.variables["XLAT"][str_t:end_t, str_i:end_i, str_j:end_j]
lat_units = "south_north"
lon = file.variables["XLONG"][str_t:end_t, str_i:end_i, str_j:end_j]
lon_units = "west_east"

ph = file.variables["PH"][str_t:end_t, :, str_i:end_i, str_j:end_j]
phb = file.variables["PHB"][str_t:end_t, :, str_i:end_i, str_j:end_j]
pb = file.variables["PB"][str_t:end_t, :, str_i:end_i, str_j:end_j]
p = file.variables["P"][str_t:end_t, :, str_i:end_i, str_j:end_j]
REFL_10CM = file.variables["REFL_10CM"][str_t:end_t, :, str_i:end_i, str_j:end_j]

# Derived variables
W = getvar(file, "wa", timeidx=-1)[str_t:end_t, :, str_i:end_i, str_j:end_j]
z = getvar(file, "z", timeidx=-1)[str_t:end_t, :, str_i:end_i, str_j:end_j]
tk = getvar(file, "T", timeidx=-1)[str_t:end_t, :, str_i:end_i, str_j:end_j] - 273.15

# Hydrometeor fields
qc = file.variables["QCLOUD"][str_t:end_t, :, str_i:end_i, str_j:end_j]
qr = file.variables["QRAIN"][str_t:end_t, :, str_i:end_i, str_j:end_j]
qi = file.variables["QICE"][str_t:end_t, :, str_i:end_i, str_j:end_j]
qs = file.variables["QSNOW"][str_t:end_t, :, str_i:end_i, str_j:end_j]
qg = file.variables["QGRAUP"][str_t:end_t, :, str_i:end_i, str_j:end_j]
EVAPPROD = file.variables["EVAPPROD"][str_t:end_t, :, str_i:end_i, str_j:end_j]
RAINPROD = file.variables["RAINPROD"][str_t:end_t, :, str_i:end_i, str_j:end_j]

# Gas species
o3 = file.variables["o3"][str_t:end_t, :, str_i:end_i, str_j:end_j]
no = file.variables["no"][str_t:end_t, :, str_i:end_i, str_j:end_j]
no2 = file.variables["no2"][str_t:end_t, :, str_i:end_i, str_j:end_j]
co = file.variables["co"][str_t:end_t, :, str_i:end_i, str_j:end_j]
hcho = file.variables["hcho"][str_t:end_t, :, str_i:end_i, str_j:end_j]
h2o2 = file.variables["h2o2"][str_t:end_t, :, str_i:end_i, str_j:end_j]
ch3ooh = file.variables["ch3ooh"][str_t:end_t, :, str_i:end_i, str_j:end_j]

# Offsets
hcho_off = file_off.variables["hcho"][str_t:end_t, :, str_i:end_i, str_j:end_j]
h2o2_off = file_off.variables["h2o2"][str_t:end_t, :, str_i:end_i, str_j:end_j]
ch3ooh_off = file_off.variables["ch3ooh"][str_t:end_t, :, str_i:end_i, str_j:end_j]

# Calculate derived variables
hcho_se = 100 * ((hcho_off - hcho) / hcho_off)
h2o2_se = 100 * ((h2o2_off - h2o2) / h2o2_off)
ch3ooh_se = 100 * ((ch3ooh_off - ch3ooh) / ch3ooh_off)
qtot = qc + qr + qi + qs + qg

# Save output
with Dataset(output_dir + output_filename, "w", format="NETCDF4") as fout:
    fout.createDimension("Time", times.shape[0])
    fout.createDimension("bottom_top", qc.shape[1])
    fout.createDimension("south_north", lat.shape[1])
    fout.createDimension("west_east", lon.shape[2])

    # Define variables and attributes as per your NCL script
    lat_var = fout.createVariable("XLAT", "f4", ("Time", "south_north", "west_east"))
    lat_var[:] = lat
    lat_var.units = lat_units

    # Add other variables similarly...
    print("Output saved to {}".format(output_filename))
