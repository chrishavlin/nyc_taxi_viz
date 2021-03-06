""" 
pl_plot_speed.py
plots results from pl_proc_speed.py

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
load_dir='../data_products/hysteresis_3mo'

# load saved files
N=tm.read_taxi_count_speed(load_dir,'N')
V=tm.read_taxi_count_speed(load_dir,'Speed')
t=tm.read_taxi_count_speed(load_dir,'t')

# get average of days
N_unAve=np.mean(N,axis=0)
SpeedAve=np.mean(V,axis=0)

# plotting
plt.figure()
plt.plot(N_unAve,SpeedAve, 'k')
plt.hold
plt.scatter(N_unAve,SpeedAve, c=t,s=50,edgecolor='none')
plt.colorbar()
for it in range(0,24,1):
    it_id=np.where(abs(t-it)==min(abs(t-it)))
    if len(it_id[0])>1:
       it_id=it_id[0][0]
    else:
       it_id=it_id[0]
    plt.scatter(N_unAve[it_id],SpeedAve[it_id],s=125,c='k')

#plt.figure()
#Ndays=len(V[:,0])
#for iday in range(0,Ndays):
#    clr=float(iday)/float(Ndays)
#    plt.scatter(N[iday,:],V[iday,:],c=(clr,0,0))
#    plt.title('day'+str(iday)+'of'+str(Ndays))
#    plt.pause(0.05)
     
     
plt.show()

