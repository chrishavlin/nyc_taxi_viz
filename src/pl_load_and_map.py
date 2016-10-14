""" 
pl_load_and_map.py
a script for loading and plotting NYC yellowcab data already processed and binned in lat/lon

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
import taxi_plotmod as tpm
import taxi_main as tm
import numpy as np
import matplotlib.pyplot as plt

# select directory and variable
data_dir='../data_products/sub_sampled'
varname='dist_mi'

# output figure names
fig_name='dist_map'

# load binned data 
VarMean,VarCount,Varx,Vary=tm.read_gridded_file(data_dir,varname)   
VarMean=VarMean * (VarCount > 5)

# plot settings
log_plot=False # take log10 of the data?
minval=0.1 # min value for color scale
maxval=10 # max value for color scale
x_dim=5 # figure x dimension [inches]
y_dim=7 # figure y dimension [inches]
ShowFig=True # display the figure?
SaveFig=True # save the figure?
figdpi=1200 # dpi resolution for save
savename=data_dir+'/'+fig_name # name of figure for save

# plot it
tpm.plt_map(VarMean,minval,maxval,Varx,Vary,log_plot,ShowFig,SaveFig,savename,x_dim,y_dim,figdpi)


# calculate and plot mean trip distance vs latitude
Mean_vs_y = np.zeros((Vary.size,1))

for iLat in range(0,len(Vary)-1):
    x_var=VarMean[iLat,:]
    x_var=x_var[np.where(x_var>0.1)]
    if len(x_var)!=0:
       Mean_vs_y[iLat] = np.mean(x_var)

fig=plt.figure()
plt.plot(Mean_vs_y,Vary)
plt.xlabel('Mean trip distance in miles')
plt.ylabel('Degree latitude')
axes = plt.gca()
axes.set_ylim([Vary.min(),Vary.max()])
axes.set_xlim([2,4.5])

print '...'
print '... close figure to continue ...'
plt.show()

