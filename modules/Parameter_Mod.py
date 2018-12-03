

pts = True  # print some informations to screen

interpolate_cn = False
create_nc_file = False
calc_doppler_spectra = False
plot_radar_results   = False
plot_compare_noise   = False
plot_for_poster      = False
plot_comparisons     = False
plot_interp2d        = False
plot_interpolation_scatter = False

# constants
chirpTable_min_height = 100.0    # (m)

# minimum and maximum height for plotting
h_min = 0.0
h_max = 15.0

# dots per inch parameter for .png
dpi_val = 200

# file extensions of different NetCDF datasets
LIMRAD_file_extension = '*.LV1.NC'
mira_file_extension   = '*mira.nc'
mmclx_file_extension  = '*.mmclx'

# path to meteorological data
meteo_path  = '/Users/willi/data/MeteoData/'
# LIMRAD_path = '/Volumes/Data_Storag/MeteoData/LIMRad94/LV0/'

LIMRAD_path = '/Users/willi/data/MeteoData/LIMRad94/noise/'
# LIMRAD_path = '/Users/willi/data/MeteoData/LIMRad94/calibrated/all/LV1'
MIRA_path   = '/Users/willi/data/MeteoData/MIRA/'

#interp_meth = 'linear'
interp_meth = 'NearestNeighbour'

interp_time_res  = 10  # in seconds
interp_range_res = 30  # in meter
