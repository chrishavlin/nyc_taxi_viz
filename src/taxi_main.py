""" 
taxi_main.py
module for loading the raw csv taxi files.

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
import time,os 
import matplotlib.pyplot as plt
from matplotlib import cm
import taxi_plotmod as tpm

"""---------
Functions
------------"""
def read_all_variables(f,there_is_a_header,VarImportList):
    """ 
     reads in the raw data from a single file
    
     input: 
       f                   file object
       there_is_a_header   logical flag 
       VarImportList       a list of strings identifying which
                           data to read in and save
               possible variables: 'pickup_time_hr','dist_mi','speed_mph','psgger','fare',
                                   'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
                                   'drop_lat','elapsed_time_min'
       
     output: 
       Vars                a 2D array, each row is a single taxi
                           pickup instance, each column is a different
       Var_list            a list of strings where the index of each 
                           entry corresponds to the column of Vars

    """

#   count number of lines
    indx=0
    for line in f:
        indx=indx+1
    if there_is_a_header:
        indx = indx-1
    Nlines = indx
#   Initizialize Variable Array and List
    Vars=np.zeros((indx,len(VarImportList)))
    Var_list=[None] * len(VarImportList)

#   Go back to start of file, loop again to read variables
    f.seek(0)
    if there_is_a_header:
        headerline=f.readline()

    indx=0

#   loop over lines, store variables
    prevprog=0
    for line in f:

        prog= round(float(indx) / float(Nlines-1) * 100)
        if prog % 5 == 0 and prog != prevprog and Nlines > 500:
           print '    ',int(prog),'% of file read ...'
           prevprog=prog

        line = line.rstrip()
        line = line.split(',')

        var_indx = 0        
        if len(line) == 19:

           if 'pickup_time_hr' in VarImportList:
              Vars[indx,var_indx]=datetime_string_to_time(line[1],'hr')
              Var_list[var_indx]='pickup_time_hr'
              var_indx=var_indx+1

           if 'date' in VarImportList:
              Vars[indx,var_indx]=line[1].split()[0] # the date string, "yyyy-mm-dd"
              Var_list[var_indx]='date'
              var_indx=var_indx+1
           
           if 'dropoff_time_hr' in VarImportList:
              Vars[indx,var_indx]=datetime_string_to_time(line[2],'hr')
              Var_list[var_indx]='dropoff_time_hr'
              var_indx=var_indx+1

           if 'dist_mi' in VarImportList:
              Vars[indx,var_indx]=float(line[4]) # distance travelled [mi]
              Var_list[var_indx]='dist_mi'
              var_indx=var_indx+1
              
           if 'elapsed_time_min' in VarImportList:
              pickup=datetime_string_to_time(line[1],'min')
              drop=datetime_string_to_time(line[2],'min')
              Vars[indx,var_indx]=drop - pickup
              Var_list[var_indx]='elapsed_time_min'
              var_indx=var_indx+1

           if 'speed_mph' in VarImportList:
              pickup=datetime_string_to_time(line[1],'min')
              drop=datetime_string_to_time(line[2],'min')
              dist=float(line[4]) # [mi]
	      if (drop-pickup) > 0:
                   speed=dist / (drop - pickup) # [mi/hr]
	      else:
	           speed=0
              Vars[indx,var_indx]=speed
              Var_list[var_indx]='speed_mph'
              var_indx=var_indx+1

           if 'pickup_lat' in VarImportList:
              Vars[indx,var_indx]=float(line[6])
              Var_list[var_indx]='pickup_lat'
              var_indx=var_indx+1

           if 'pickup_lon' in VarImportList:
              Vars[indx,var_indx]=float(line[5])
              Var_list[var_indx]='pickup_lon'
              var_indx=var_indx+1

           if 'drop_lat' in VarImportList:
              Vars[indx,var_indx]=float(line[10])
              Var_list[var_indx]='drop_lat'
              var_indx=var_indx+1

           if 'drop_lon' in VarImportList:
              Vars[indx,var_indx]=float(line[9])
              Var_list[var_indx]='drop_lon'
              var_indx=var_indx+1

           if 'psgger' in VarImportList:
              Vars[indx,var_indx]=float(line[3])
              Var_list[var_indx]='pssger'
              var_indx=var_indx+1

           if 'fare' in VarImportList:
              Vars[indx,var_indx]=float(line[12])
              Var_list[var_indx]='fare'
              var_indx=var_indx+1

           if 'tips' in VarImportList:
              Vars[indx,var_indx]=float(line[15])
              Var_list[var_indx]='tips'
              var_indx=var_indx+1

           if 'payment_type' in VarImportList:
              Vars[indx,var_indx]=float(line[11])
              Var_list[var_indx]='payment_type'
              var_indx=var_indx+1

           indx=indx+1

    return Vars,Var_list

def datetime_string_to_time(dt_string,time_units):
    """ converts datetime string to time in units of time_units 
        dt_string should be in datetime format: "yyyy-mm-dd hh:mm:ss" 
                                                "2016-04-18 18:31:43"
    
    """
    t_string=dt_string.split()[1] # remove the space, take the time string 
    t_hms=t_string.split(':') # split into hr, min, sec 

    # unit conversion factors depending on time_units:
    if time_units == 'hr':
       a = [1, 1/60, 1/3600] 
    elif time_units == 'min':
       a = [60, 1, 1/60] 
    elif time_units == 'sec':
       a = [3600, 60, 1] 

    time_flt=float(t_hms[0])*a[0]+float(t_hms[1])*a[1]+float(t_hms[2])*a[2]

    return time_flt


def read_taxi_files(dir_base,Vars_To_Import):
    """ loops over all taxi files in a directory, stores them in memory
       
        input:
           dir_base        the directory to look for .csv taxi files
           Vars_to_Import  a list of strings identifying which data to read in and save

               possible variables: 'pickup_time_hr','dist_mi','speed_mph','psgger','fare',
                                   'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
                                   'drop_lat','elapsed_time_min'
       output: 
        

       VarBig          a 2D array, each row is a single taxi pickup instance, each column
                       is a different variable. Data aggregated from all files in directory.
       Var_list         a list of strings where the index of each entry corresponds to the 
                        column of Vars
    
    """

    N_files=len(os.listdir(dir_base)) # number of files in directory
    ifile = 1 # file counter
    Elapsed_tot=0 # time counter
    
    for fn in os.listdir(dir_base):             # loop over directory contents 
         if os.path.isfile(dir_base+fn):        # is the current path obect a file?
            flnm=dir_base + fn                  # construct the file name
            print 'Reading File ', ifile,' of ', N_files
    
            start = time.clock()                # start timer
            fle = open(flnm, 'r')               # open the file for reading
    
          # distribute current file to lat/lon bins:

            VarChunk,Var_list=read_all_variables(fle,True,Vars_To_Import) 
    
            if ifile == 1:
    	       VarBig = VarChunk
    	    else: 
               VarBig = np.vstack((VarBig,VarChunk))
    
            elapsed=(time.clock()-start)        # elapsed time
            Elapsed_tot=Elapsed_tot+elapsed     # cumulative elapsed
            MeanElapsed=Elapsed_tot/ifile       # mean time per file
            Fls_left=N_files-(ifile)            # files remaining 
            time_left=Fls_left*MeanElapsed/60   # estimated time left
    
            print '    aggregation took %.1f sec' % elapsed
            print '    estimated time remaning: %.1f min' % time_left
    
            fle.close()                         # close current file
            ifile = ifile+1                     # increment file counter

    return VarBig,Var_list

def write_gridded_file(write_dir,Var,VarCount,x,y,Varname):   
    """ writes out the spatially binned data """

    if not os.path.exists(write_dir):
       os.makedirs(write_dir)

    f_base=write_dir+'/'+Varname
    np.savetxt(f_base +'.txt', Var, delimiter=',')
    np.savetxt(f_base +'_Count.txt', VarCount, delimiter=',')
    np.savetxt(f_base+'_x.txt', x, delimiter=',')
    np.savetxt(f_base+'_y.txt', y, delimiter=',')

def read_gridded_file(read_dir,Varname):   
    """ reads in the spatially binned data """
    f_base=read_dir+'/'+Varname
    Var=np.loadtxt(f_base +'.txt',delimiter=',')
    VarCount=np.loadtxt(f_base +'_Count.txt',delimiter=',')
    x=np.loadtxt(f_base+'_x.txt',delimiter=',')
    y=np.loadtxt(f_base+'_y.txt',delimiter=',')

    return Var,VarCount,x,y

""" END OF FUNCTIONS """

if __name__ == '__main__':
    """ a basic example of reading, processing and plotting some taxi files """

    # the directory with the data
    dir_base='../data_sub_sampled/'

    # choose which variables to import
    # possible variables: 'pickup_time_hr','dist_mi','speed_mph','psgger','fare',
    #          'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
    #          'drop_lat','elapsed_time_min'

    Vars_To_Import=['dist_mi','pickup_lon','pickup_lat']

    # read in all the data! 
    VarBig,Var_list=read_taxi_files(dir_base,Vars_To_Import)
    
    # now bin the point data!
    DistCount,DistMean,Distx,Disty=tpm.map_proc(VarBig,Var_list,'dist_mi',0.1,60,'True',600,700)
    write_gridded_file('../data_products/',DistMean,DistCount,Distx,Disty,'dist_mi')   

    tpm.plt_map(DistCount,1,1000,Distx,Disty,True)

