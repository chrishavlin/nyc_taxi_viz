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

nt = int(24*4)
time_grid = np.linspace(0,24,nt)

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
write_dir='./data_timeN__' + str(nt)
try: 
    os.makedirs(write_dir)
except OSError:
    if not os.path.isdir(write_dir):
       raise


"""----------------- count_taxi_file ----------------------- """
def count_taxi_file(f,nt,time_grid,there_is_a_header):
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
    TaxiCount=np.zeros((nt-1,1)) # defined on cell centers

#   count number of lines, define the latitude, longitude arrays
    indx=0
    for line in f:
        indx=indx+1
    
    if there_is_a_header:
        indx = indx-1

    BigTime=np.zeros((indx,1))

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
	   pickup=pickup.split() # remove the space
	   pickup=pickup[1] # take time only (date is pikcup[0])
	   t_hms=pickup.split(':')
	   pickup_time=float(t_hms[0])+float(t_hms[1])/60+float(t_hms[2])/3600

           BigTime[indx,0]=pickup_time
           indx=indx+1
    
#   count the number of points within each time bin    
    for it in range(nt-1):
            tmin=time_grid[it] 
            tmax=time_grid[it+1] 
            
	    BigTime2=BigTime[BigTime[:,0]>tmin]
	    BigTime2=BigTime2[BigTime2[:,0]<tmax]

            TaxiCount[it]=len(BigTime2)

    return TaxiCount
"""-------------- END count_taxi_file ----------------------- """

"""----------------------------------------------------------
Processing: 
Loops over data files, counts number of taxi pickups within
each lat/lon bin
-------------------------------------------------------------"""
TaxiCount=np.zeros((nt-1,1)) # zero array to add counts to
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
        TaxiCount2=count_taxi_file(fle,nt,time_grid,True) 

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


time_t = (time_grid[0:nt-1] + time_grid[1:nt])/2
plt.plot(time_t,TaxiCount)
plt.show()
"""-------------------------- end of taxi_map_plot.py ----------------------"""
