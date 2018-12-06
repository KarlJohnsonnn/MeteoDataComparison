# Descibtion

This is the version 0.1 of the Meteorological-Data-Comparision Package. It's perpous is to investigate RPG 94GHz FMCW Radar 'LIMRAD94' and Metek 35GHz Pulse Radar 'MIRA35' radar data files. The user is able generates quicklooks for quick investigation of radar moments (Ze, mdv, sw, ldr), looking at specific Doppler spectras, plotting them, convert LIMRAD94 data to Cloudnet processing format, etc.


# Installation

  1.  make sure to install Python 3.6 or a later verison (use anaconda3 for your discribution)
  
  2.  the following packages are nessessarry ( type:  pip install [packagename] or conda install [packagename])
        - numpy
        - matplotlib
        - netCDF4
        - numba
        - scipy
        - ...
        
  3. Define global paths via Parameter_Mod.py file in the subfolder modules/. 
     Make a copy of Parameter_Mod.py_untouched and save it as Parameter_Mod.py. 
     
    $ cp Parameter_Mod.py_untouched Parameter_Mod.py
    
     Then specify your local paths:
     
    - meteo_path  = '/Users/willi/data/MeteoData/'              # path to radar data home path
    - LIMRAD_path = '/Users/willi/data/MeteoData/LIMRAD94/'     # path to LIMRAD94 NetCDF files
    - MIRA_path   = '/Users/willi/data/MeteoData/MIRA/'         # path to MIRA NetCDF files
  
