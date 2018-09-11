

pts = True # print to screen
dbg = False # debug ouput flag


# constants
chirpTable_min_height = 0.1 # (km)

# dots per inch parameter for .png
dpi_val = 200

# file extensions of different NetCDF datasets
LIMRad_file_extension = '*.LV1.NC'
mira_file_extension   = '*mira.nc'
mmclx_file_extension  = '*.mmclx'


plot_interp2d              = False
plot_RectBivariateSpline   = False           #interp2 testen
plot_radar_results         = False
plot_comparisons           = False
plot_interpolation_scatter = False
plot_compare_mira_mmclx    = False


# path to meteorological data
meteo_path = '/Users/willi/MeteoData/'

interp_meth = 'NearestNeighbour'
