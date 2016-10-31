""" 
pl_proc_speed.py
processes the raw data to extract taxi speed and number of taxis on road

Copyright (C) 2016  Chris Havlin, <https://chrishavlin.wordpress.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. 

The database is NOT distributed with the code here. 

Data source:
     NYC Taxi & Limousine Commision, TLC Trip Record Data
     <http://www.nyc.gov/html/tlc/html/about/trip_record_data.shtml>
"""

"""--------------
Import libraries:
-----------------"""
import numpy as np
import taxi_plotmod as tpm
import taxi_main as tm
import matplotlib.pyplot as plt
import datetime as dt


# the directory with the data
#dir_base='../data_sub_sampled/'
dir_base='../full_csv_files/'

# choose which variables to import
# possible variables: 'pickup_time_hr','dist_mi','speed_mph','psgger','fare',
#          'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
#          'drop_lat','elapsed_time_min'

# import data
Vars_To_Import=['elapsed_time_min','pickup_time_hr','speed_mph']
VarBig,Var_list,Date=tm.read_taxi_files(dir_base,Vars_To_Import)

# cull data 
print 'shape, before culling:',VarBig.shape
# speed limit
Speed=VarBig[:,Var_list.index('speed_mph')]
ids=np.where((VarBig[:,Var_list.index('speed_mph')]<80) &
              (VarBig[:,Var_list.index('speed_mph')]>0))
VarBig=VarBig[ids[0],:]
Date=Date[ids[0]]

# elapsed time limit
Ela=VarBig[:,Var_list.index('elapsed_time_min')]
idA=np.where((Ela<120) & (Ela > 0))
VarBig=VarBig[idA[0],:]
Date=Date[idA[0]]
print 'shape, after culling:',VarBig.shape

# loop over dates
date_start=min(Date)
date_end=max(Date)
date_i = date_start
one_day_delta = np.timedelta64(1, 'D')

# calculate first date
ids=np.where(Date==date_start)
VarDate=VarBig[ids[0],:]#tpm.select_data_by_date(VarBig,Var_list,Date,date_i.year,date_i.month,date_i.day)
N_unAve,SpeedAve,t_unAve = tpm.find_N_unique_vs_t(VarDate,Var_list)
date_count=1
date_i=date_i+one_day_delta

print date_start,date_end,type(date_i)

while date_i <= date_end:
      print 'processing',date_i

      # process single day
      ids=np.where(Date==date_i)
      VarDate=VarBig[ids[0],:]#tpm.select_data_by_date(VarBig,Var_list,Date,date_i.year,date_i.month,date_i.day)
#      VarDate=tpm.select_data_by_date(VarBig,Var_list,date_i.year,date_i.month,date_i.day)
      print 'Shapez:',VarBig.shape,VarDate.shape
      N_un,Speed_un,t_un = tpm.find_N_unique_vs_t(VarDate,Var_list)
#
      # keep track of sum
      N_unAve=N_unAve+N_un
      SpeedAve=SpeedAve+Speed_un
#
      # keep only unprocessed dates
      ids=np.where(Date>date_i)
      VarBig=VarBig[ids[0],:]
      Date=Date[ids[0]]

      # increment counters
      date_count=date_count+1
      date_i = date_i + one_day_delta
#
# get average of days
N_unAve=N_unAve/float(date_count)
SpeedAve=SpeedAve/float(date_count)

## do it all in one chunk
#N_unAve2,SpeedAve2,t_unAve2 = tpm.find_N_unique_vs_t(VarBig,Var_list)
#

plt.figure()
plt.plot(N_unAve,SpeedAve, 'k')
plt.scatter(N_unAve,SpeedAve, c=t_unAve)
plt.colorbar()
plt.show()

