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


# the directory with the data
dir_base='../data_sub_sampled/'

# choose which variables to import
# possible variables: 'pickup_time_hr','dist_mi','speed_mph','psgger','fare',
#          'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
#          'drop_lat','elapsed_time_min'

Vars_To_Import=['date','speed_mph','pickup_lon','pickup_lat']

# read in all the data! 
VarBig,Var_list=tm.read_taxi_files(dir_base,Vars_To_Import)
    

print VarBig[0:10,0]    
## now bin the point data and save the result!
#x_bins=600 # number of bins in x
#y_bins=700 # number of bins in y

#DistCount,DistMean,Distx,Disty=tpm.map_proc(VarBig,Var_list,'dist_mi',0.1,60,True,x_bins,y_bins)
#tm.write_gridded_file('../data_products/sub_sampled/',DistMean,DistCount,Distx,Disty,'dist_mi')   

#FareCount,FareMean,Farex,Farey=tpm.map_proc(VarBig,Var_list,'fare',0.1,60,True,x_bins,y_bins)
#tm.write_gridded_file('../data_products/sub_sampled/',FareMean,FareCount,Farex,Farey,'fare')   
