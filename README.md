# Descibtion

This is the version 0.1 of the Meteorological-Data-Comparision Package. It's perpous is to investigate RPG 94GHz FMCW Radar 'LIMRAD94' and Metek 35GHz Pulse Radar 'MIRA35' radar data files. The user is able generates quicklooks for quick investigation of radar moments (Ze, mdv, sw, ldr), looking at specific Doppler spectras, plotting them, convert LIMRAD94 data to Cloudnet processing format, etc.


# Installation

  1.  make sure to install Python 3.6 or a later verison, e.g. go to https://www.anaconda.com/download/ and pick a download depending on your operating system (this may take a while)
      
  
  2.  the following packages are nessessarry, type: conda install [packagename]
     
     $ conda install -c anaconda numpy
     $ conda install -c conda-forge matplotlib
     $ conda install -c conda-forge netcdf4
     $ conda install -c anaconda numba 
     $ conda install -c anaconda scipy 
        
  3. Make a copy of Parameter_Mod.py_untouched and save it as Parameter_Mod.py. 
     Then define global paths in Parameter_Mod. 
     
    $ cp Parameter_Mod.py_untouched Parameter_Mod.py
    
     Then specify your local paths, e.g.:
     
    - meteo_path  = '[user]/data/MeteoData/'              # path where output is stored, e.g.: png, log, txt
    - LIMRAD_path = '[user]/data/MeteoData/LIMRAD94/'     # main path to LIMRAD94 NetCDF files
    - MIRA_path   = '[user]/data/MeteoData/MIRA/'         # main path to MIRA NetCDF files
   
   The folder structure is as follows:
    
      MIRA_parth/mmclx/[files].mmclx
         --"--  /calibrated/[files].mira
         --"--  /spectra/[files].nc4
                 
      /LIMRAD_path/calibrated/[momentfiles].LV1.NC   
                --"--        /[spectrafiles].LV0.NC           
          
  
