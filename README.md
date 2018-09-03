# LIMRAD94-MIRA35-comparison.py
This is version 0.1 of "LIMRAD94-MIRA35-comparison" and has the perpous of the quick investigation of differnces between the RPG 94GHz FMCW Radar 'LIMRad94' and the Metek 35GHz Pulse Radar 'MIRA35'.

The analysis is done by comparing three main radar moments: 
  - equivalent radar reflectivity factor (Ze)
  - mean Doppler velocity (mdv)
  - spectral width (sw)

The input routine is able to handle ".LV1.NC" NetCDF data (LIMRad94) as well as ".mira" and ".mmclx" NetCDF data (MIRA35).


Installation:

  1.  make sure to install Python 3.6 or a later verison
  
  2.  the following packages are nessessarry ( type:  pip install [packagename] )
        - numpy
        - matplotlib
        - netCDF4
  
Execution:
  
  The user has to execute the code with five additional arguments in order to specify the exact time and location of an event:
    
    - minheight ... (float) minimum height of the desired frame, e.g.: 8.5
    - maxheight ... (float) maximum height of the desired frame, e.g.: 12.0
    - date ........ (integer) date of the dataset format: YYMMDD, e.g.: 180728
    - from ........ (integer) starting time of the desired date in UTC, format: HHMM, e.g.: 0710
    - to .......... (integer) end time of the desired date in UTC, format: HHMM, e.g.: 1240

  $ python3 LIMRAD94-MIRA35-comparison.py minheight maxheight date from to


Program Settings:

  The Python code contains several logical variables, which allow the execution of different tasks.
  These parameters have to be set by the user before execution, in order to output the desired information and graphics.
  If a parameter is equal to True, a figure containing different subplots is created and stored as .png.

  possible scenario:

    - plot_RectBivariateSpline   = False    # interpolates 2D radar reflectivity data (experimental)
    - plot_radar_results         = True     # plotting the radar moments (Ze,mdv,sw) of LIMRad94 and MIRA35 NetCDF data
    - plot_comparisons           = True     # computes hight and time-averaged data of LIMRad94 and MIRA35 and plots the results
    - plot_interpolation_scatter = True     # interpolating mean-height onto a uniformly spaced grid, calulation of means and correlation coefficient, scatter plot
    - plot_compare_mira_mmclx    = True     # comparision of different MIRA35 modes (considering all targets, only hydrometeors, ...)

    - pts = True                          # print some information to screen (progress, file names of .png output)
    - dbg = False                         # if True some more information is printed to screen (array dimensions, ...)

    - LIMRad_file_extension = '*.LV1.NC'  # this is a fixed parameter, do not alter
    - mira_file_extension   = '*.mmclx'   # either:  '*mira.nc'   processed data, or:  '*.mmclx'   processed and unprocessed data

    - chirpTable_min_height = 0.1         # minimum height of the first chirp sequence of LIMRad94 in [km], will be automated in the future
  


