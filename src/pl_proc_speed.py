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

write_dir='../data_products/hysteresis_test'

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

# plot 2-D histogram
time_bin_edge=np.linspace(0,24,24*4)
tpm.plt_two_d_histogram('speed_mph',0,80,time_bin_edge,VarBig,Var_list)

# loop over dates
date_start=min(Date)
date_end=max(Date)
date_i = date_start
one_day_delta = np.timedelta64(1, 'D')

t_unAve=np.linspace(0,24,24*10)

# find total elapsed days, initialize N and Speed arrays
Ndates=date_end.astype(dt.datetime)-date_start.astype(dt.datetime) # (convert from np datetime64 to datetime)
N_un=np.zeros((Ndates.days+1,len(t_unAve)))
Speed_un=np.zeros((Ndates.days+1,len(t_unAve)))

# now loop through days to find daily hysterises 
date_count=0
while date_i <= date_end:
      print 'processing',date_i

      # pull out and process single day
      ids=np.where(Date==date_i)
      VarDate=VarBig[ids[0],:]
      N_un[date_count,:],Speed_un[date_count,:],t_un = tpm.find_N_unique_vs_t(VarDate,Var_list,t_unAve)

      # keep only unprocessed dates
      ids=np.where(Date>date_i)
      VarBig=VarBig[ids[0],:]
      Date=Date[ids[0]]

      # increment counters
      date_count=date_count+1
      date_i = date_i + one_day_delta

# write out the Speed and N and t arrays
tm.write_taxi_count_speed(write_dir,N_un,'N',Speed_un,'Speed',t_unAve,'t')

# test the load
N=tm.read_taxi_count_speed(write_dir,'N')
V=tm.read_taxi_count_speed(write_dir,'Speed')
t=tm.read_taxi_count_speed(write_dir,'t')

# get average of days
N_unAve=np.mean(N,axis=0)
SpeedAve=np.mean(V,axis=0)

# plotting
plt.figure()
plt.plot(N_unAve,SpeedAve, 'k')
for it in range(0,24,1):
    it_id=np.where(abs(t-it)==min(abs(t-it)))
    plt.scatter(N_unAve[it_id[0]],SpeedAve[it_id[0]],s=125,c='k')
plt.scatter(N_unAve,SpeedAve, c=t,s=50,edgecolor='none')
plt.colorbar()

plt.show()

