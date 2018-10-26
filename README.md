# Playground.py
This is the main file version 0.1 of "LIMRAD94-MIRA35-comparison" and has the perpous of showing the user how to input/ouput LIMRAD94 and MIRA NetCDF radar files. It generates quicklooks for the investigation of radar moments (Ze, mdv, sw, ldr) messured by the RPG 94GHz FMCW Radar 'LIMRAD94' and the Metek 35GHz Pulse Radar 'MIRA35'.

The analysis is done by comparing three main radar moments: 
  - equivalent radar reflectivity factor (Ze)
  - mean Doppler velocity (mdv)
  - spectral width (sw)
  - linear depolarization ratio (ldr)
  
Calculating the noise floor for uncompressed LIMRAD94 .LV0.NC and compaEna = 0 files.

The input routine is able to handle ".LV0.NC" and ".LV1.NC" NetCDF data (LIMRad94) as well as ".mira" and ".mmclx" NetCDF data (MIRA35). The Playground.py should serve the user how to call the input, output and plotting routines. 


Installation:

  1.  make sure to install Python 3.6 or a later verison
  
  2.  the following packages are nessessarry ( type:  pip install [packagename] )
        - numpy
        - matplotlib
        - netCDF4
        - numba
        - scipy
  
Execution:
  
  The user has to execute the code with five additional arguments in order to specify the exact time and location of an event:
    
    - date ........ (integer) date of the dataset format: YYMMDD, e.g.: 180728                    (default: 180729)
    - from ........ (integer) starting time of the desired date in UTC, format: HHMM, e.g.: 0710  (default:   0000)
    - to .......... (integer) end time of the desired date in UTC, format: HHMM, e.g.: 1240       (default:   2359)
    - minheight ... (float) minimum height of the desired frame, e.g.: 8.5                        (default:    0.0)
    - maxheight ... (float) maximum height of the desired frame, e.g.: 12.0                       (default:   12.0)


# Calling the routine with specific parameters

    $ python Playground.py minheight maxheight date from to


# Available .png Outputs

  The Python code contains several logical variables, which allow the execution of different tasks.
  These parameters have to be set by the user before execution, in order to output the desired information and graphics.
  If a parameter is equal to True, a figure containing different subplots is created and stored as .png.

  possible scenario:

    - plot_radar_results         = True     # plotting the radar moments (Ze,mdv,sw) of LIMRad94 and MIRA35 NetCDF data
    - plot_comparisons           = True     # computes hight and time-averaged data of LIMRad94 and MIRA35 and plots the results
    - plot_interpolation_scatter = True     # interpolating mean-height onto a uniformly spaced grid, calulation of means and correlation coefficient, scatter plot
    - plot_for_poster            = False    # plots Ze and mdv with larger fontsize
    - plot_interp2d              = False    # plots interpolation of Ze from LIMRAD94<->MIRA grid, calculates correlation
    - #plot_compare_mira_mmclx   = False    # disabled 
    - plot_doppler_spectra       = False    # plots Doppler spectra from Lv0 files with noise floor estimation

# Parameter_Mod.py
  Contains some global variabels and flags (default values below):
  
    - pts = True                          # print some information to screen (progress, file names of .png output)
    
    - chirpTable_min_height = 0.1         # minimum height of the first chirp sequence of LIMRad94 in [km], will be automated in the future
       
    - LIMRAD_lv0_fext = '*.LV0.NC'        # LIMRAD94 Lv0 datafile (moments)
    - LIMRAD_lv1_fext = '*.LV1.NC'        # LIMRAD94 Lv1 datafile (moments)
    - mira_fext  = '*mira.nc'             # processed data of MIRA35 GHz Radar
    - mmclx_fext = '*.mmclx'              # contains processed and unprocessed of MIRA35 GHz Radar
  
    - meteo_path  = '/Users/willi/data/MeteoData/'              # path to radar data home path
    - LIMRAD_path = '/Users/willi/data/MeteoData/LIMRAD94/'     # path to LIMRAD94 NetCDF files
    - MIRA_path   = '/Users/willi/data/MeteoData/MIRA/'         # path to MIRA NetCDF files

    - interp_meth = 'NearestNeighbour'    # interpolation method for regrid ('NearestNeighbour', 'linear')
    - interp_time_res = 10                # time resolution in seconds
    - interp_range_res = 30               # range resolution in meter
    - interpolate_cn = True               # interpolates LIMRAd94 onto a user definded grid
    - create_nc_file = True               # creates concatinated (multiple successive files) with original and interpolated resolution 
    - dpi_val = 400                       # dots per inch parameter for output .png
