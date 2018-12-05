
pts = True  # print some informations to screen

# constants
chirpTable_min_height = 100.0    # (m)

# minimum and maximum height for plotting
h_min = 0.0
h_max = 15.0

# dots per inch parameter for .png
dpi_val = 100

# file extensions of different NetCDF datasets
LIMRAD_file_extension = '*.LV*.NC'
mira_file_extension   = '*mira.nc'
mmclx_file_extension  = '*.mmclx'

# path to meteorological data
meteo_path  = '/Users/willi/data/MeteoData/'                 # path to folder where output (png, txt) is stored
LIMRAD_path = '/Users/willi/data/MeteoData/LIMRad94/noise/'  # path to folder where LIMRAD94 LV0.NC and LV1.NC is stored
MIRA_path   = '/Users/willi/data/MeteoData/MIRA/'            # path to folder where MIRA .mira and .mmclx is stored

#interp_meth = 'linear'
interp_meth = 'NearestNeighbour'

interp_time_res  = 10  # in seconds
interp_range_res = 30  # in meter
