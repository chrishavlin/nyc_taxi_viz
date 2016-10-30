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


# the directory with the data
#dir_base='../data_sub_sampled/'
dir_base='../full_csv_files/'

# choose which variables to import
# possible variables: 'pickup_time_hr','dist_mi','speed_mph','psgger','fare',
#          'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
#          'drop_lat','elapsed_time_min'

#Vars_To_Import=['date','pickup_time_hr','elapsed_time_min','speed_mph']
Vars_To_Import=['date','elapsed_time_min','pickup_time_hr','speed_mph']
VarBig,Var_list=tm.read_taxi_files(dir_base,Vars_To_Import)
Pickups=VarBig[:,Var_list.index('pickup_time_hr')]
Speed=VarBig[:,Var_list.index('speed_mph')]

#plt.subplot(1,2,1)
#plt.hist(Pickups,bins=24,label='full set',histtype='step')
#ids=np.where((VarBig[:,Var_list.index('speed_mph')]<80) &
#              (VarBig[:,Var_list.index('speed_mph')]>0))
#plt.hist(Pickups[ids[0]],bins=24,label='limited set (0 < speed < 80)',histtype='step')
#
#plt.subplot(1,2,2)
#plt.hist(Speed,bins=20,label='full set',histtype='step',range=(-10,80))
#plt.hist(Speed[ids[0]],bins=20,label='limited set',histtype='step')
#plt.show()
    

## Ntaxi vs time of day single day
#
ids=np.where((VarBig[:,Var_list.index('speed_mph')]<80) & 
              (VarBig[:,Var_list.index('speed_mph')]>0))
VarBig=VarBig[ids[0],:]

print VarBig[:,Var_list.index('speed_mph')].max()
print VarBig[:,Var_list.index('speed_mph')].min()
#
date_start,date_end=tpm.min_max_date(VarBig,Var_list)
date_select=date_start.split('-')
print date_start,date_end
#
date_select[2]=str(int(date_select[2])+1)
VarDate=tpm.select_data_by_date(VarBig,Var_list,date_select[0],date_select[1],date_select[2])
N_unAve,SpeedAve,t_unAve = tpm.find_N_unique_vs_t(VarDate,Var_list)


idays=1.0

for iday in range(1,30):
    print iday
    date_select[2]=str(int(date_select[2])+1)
    VarDate=tpm.select_data_by_date(VarBig,Var_list,date_select[0],date_select[1],date_select[2])
    N_un,Speed,t_un = tpm.find_N_unique_vs_t(VarDate,Var_list)

    N_unAve=N_unAve+N_un
    SpeedAve=Speed+SpeedAve
    idays=idays+1.0 

    plt.subplot(2,1,1)
    plt.plot(t_un,N_un)
    plt.xlim((21,24.25))
    plt.subplot(2,1,2)
    plt.plot(t_un,Speed)
    plt.xlim((21,24.25))
    plt.pause(0.5)

    # remove rows already processed
    ids=np.where(VarBig[:,Var_list.index('date_dd')]>iday)
    VarBig=VarBig[ids[0],:]
    print VarBig.shape

SpeedAve=SpeedAve/idays
N_unAve=N_unAve/idays

plt.subplot(2,1,1)
plt.plot(t_un,N_unAve,'k',linewidth=2.0)
plt.subplot(2,1,2)
plt.plot(t_un,SpeedAve,'k',linewidth=2.0)
plt.show()
# single day, don't need that extra info
#Var_list_r=list(Var_list)
#VarDate=np.delete(VarDate,0,Var_list_r.index('date_yr'))
#Var_list_r.remove('date_yr')
#VarDate=np.delete(VarDate,0,Var_list_r.index('date_mm'))
#Var_list_r.remove('date_mm')
#VarDate=np.delete(VarDate,0,Var_list_r.index('date_dd'))
#Var_list_r.remove('date_dd')
#
#print Var_list
#print Var_list_r


plt.scatter(N_unAve,SpeedAve, c=t_un)
plt.colorbar()
plt.show()

