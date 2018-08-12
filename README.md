# LIMRAD94-MIRA35-comparison.py
This program has the perpous of the fast investigation of differnces between the RPG 94GHz FMCW Radar 'LIMRad94' and the Metek 35GHz Pulse Radar 'MIRA35'.

The analysis is done by comparing three main radar moments: 
  - equivalent radar reflectivity factor (Ze)
  - mean Doppler velocity (mdv)
  - spectral width (sw)

The input routine is able to handle ".LV1.NC" NetCDF data (LIMRad94) as well as ".mira" and ".mmclx" NetCDF data (MIRA35).


Installation:

  1.  make sure to install Python 3.6 or a later verison
  
  2.  the following packages are nessessarry ( type:  pip3 install [packagename] )
        - numpy
        - matplotlib
        - netCDF4
  
Execution:
  
  The user has to execute the code with five additional arguments:
    
    - minheight ... (float) minimum height of the desired frame, e.g.: 8.5
    - maxheight ... (float) maximum height of the desired frame, e.g.: 12.0
    - date ........ (integer) date of the dataset format: YYMMDD, e.g.: 180728
    - from ........ (integer) starting time of the desired date in UTC, format: HHMM, e.g.: 0710
    - to .......... (integer) end time of the desired date in UTC, format: HHMM, e.g.: 1240

  $ python3 LIMRAD94-MIRA35-comparison.py minheight maxheight date from to
  

