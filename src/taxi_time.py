"""---------------------------------------------------------------------------
taxi_map_process.py

Processes the csv taxi files. 

(c) Chris Havlin
Open source license?
----------------------------------------------------------------------------"""

"""--------------
Import libraries:
-----------------"""
import numpy as np
import time,os # for 
#import io
import matplotlib.pyplot as plt
#import os # for reading/writing and directory operations

"""-----------------------------------------------------
Initialization: 
define lat/long bin grid, set read and write directories
--------------------------------------------------------"""

# extent of lat/long grid

# full manhattan (and a little brooklyn) extent
bottom_left=[40.697089, -74.027397]
top_left=[40.823067, -74.027397]
bottom_right=[40.697089, -73.914240]
top_right=[40.823067,-73.914240]

nt = int(24*2)
time_grid = np.linspace(0,24,nt)
day_array = np.linspace(1,31,31)
ndays=len(day_array)

# wall street extent
#bottom_left=[40.701035, -74.022382]
#top_left=[40.715075, -74.022382]
#bottom_right=[40.701035, -74.000207]
#top_right=[40.715075, -74.000207]

# columbus circle extent
#bottom_left=[40.765505, -73.985891]
#top_left=[40.770369, -73.985891]
#bottom_right=[40.765505,-73.978744]
#top_right=[40.770369,-73.978744]

# the directory with the data
#dir_base='./text_data/sub_sampling/'
#dir_base='./text_data/single_file/'
dir_base='./text_data/'

# make the write directory 
write_dir='./data_time_nt_' + str(nt) + '_daily'
try: 
    os.makedirs(write_dir)
except OSError:
    if not os.path.isdir(write_dir):
       raise


"""----------------- count_taxi_file ----------------------- """
def count_taxi_file(f,nt,time_grid,day_array,there_is_a_header):
    """ 
     define count_taxi_file function to count number of records
     within each lat/lon bin for a single file (to call within
     a loop over files)
    
     input: 
       f      file object
       nx     number of longitude points
       ny     number of latitutde points
       x      the x-grid array (1-D array)
       y      the y-grid array (1-D array)
       
     output: 
       TaxiCount   2D array with number of taxi pickups in each
                   lat/lon bin
    """
#   create an empty TaxiCount array 
    ndays=len(day_array)
    TaxiCount=np.zeros((nt-1,ndays)) # defined on cell centers

#   count number of lines, define the latitude, longitude arrays
    indx=0
    for line in f:
        indx=indx+1
    
    if there_is_a_header:
        indx = indx-1

    BigTime=np.zeros((indx,1))
    BigDay=np.zeros((indx,1))

#   go back to start of file, loop again to read in lat/lon    
    f.seek(0)

    if there_is_a_header:
        headerline=f.readline()

    indx=0
    for line in f:
    #   read in each line separately, extract longitude, latitutde
        line = line.rstrip()
        line = line.split(',')
        if len(line) == 19:
	   pickup=line[1] # datetime string
	   pickup=pickup.split() # split by spaces
	   
	   pickdate=pickup[0] # the date
	   pickup=pickup[1] # the time

           pickdate=int(pickdate.split('-')[2]) # day of month

	   t_hms=pickup.split(':') # split to hour, minute, second
	   pickup_time=float(t_hms[0])+float(t_hms[1])/60+float(t_hms[2])/3600 # time in hrs

           BigTime[indx,0]=pickup_time # stored
           BigDay[indx,0]=pickdate # stored
           indx=indx+1
    
#   count the number of points within each time bin    
    for it in range(nt-1):
        for iday in range(ndays):

            tmin=time_grid[it] 
            tmax=time_grid[it+1] 
	    current_day=day_array[iday]
            
	    BigDay2=BigDay[BigTime[:,0]>tmin]
	    BigTime2=BigTime[BigTime[:,0]>tmin]

	    BigDay2=BigDay2[BigTime2[:,0]<tmax]
	    BigTime2=BigTime2[BigTime2[:,0]<tmax]

            BigTime2=BigTime2[BigDay2[:,0]==current_day]
	    BigDay2=BigDay2[BigDay2[:,0]==current_day]

            TaxiCount[it,iday]=len(BigTime2)

    return TaxiCount
"""-------------- END count_taxi_file ----------------------- """

"""----------------------------------------------------------
Processing: 
Loops over data files, counts number of taxi pickups within
each lat/lon bin
-------------------------------------------------------------"""
TaxiCount=np.zeros((nt-1,ndays)) # zero array to add counts to
N_files=len(os.listdir(dir_base)) # number of files in directory
ifile = 1 # file counter
Elapsed_tot=0 # time counter

for fn in os.listdir(dir_base):             # loop over directory contents 
     if os.path.isfile(dir_base+fn):        # is the current path obect a file?
        flnm=dir_base + fn                  # construct the file name
        print 'Reading File ', ifile,' of ', N_files

        fle = open(flnm, 'r')               # open the file for reading
        start = time.clock()                # start timer

      # distribute current file to lat/lon bins:
        TaxiCount2=count_taxi_file(fle,nt,time_grid,day_array,True) 

        elapsed=(time.clock()-start)        # elapsed time
        Elapsed_tot=Elapsed_tot+elapsed     # cumulative elapsed
        MeanElapsed=Elapsed_tot/ifile       # mean time per file
        Fls_left=N_files-(ifile)            # files remaining 
        time_left=Fls_left*MeanElapsed/60   # estimated time left

        print '    aggregation took %.1f sec' % elapsed
        print '    estimated time remaning: %.1f min' % time_left

        fle.close()                         # close current file
        TaxiCount=TaxiCount+TaxiCount2      # aggregate total taxi count
        ifile = ifile+1                     # increment file counter

"""-----------------------------------------------
Save TaxiCount array and the lat/long bin arrays:
--------------------------------------------------"""

np.savetxt(write_dir+'/TaxiCount.txt', TaxiCount, delimiter=',')
np.savetxt(write_dir+'/Time.txt', time_grid, delimiter=',')
np.savetxt(write_dir+'/Days.txt', day_array, delimiter=',')


"""-------------------------- end of taxi_time.py ----------------------"""
