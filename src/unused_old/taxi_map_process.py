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
#import matplotlib.pyplot as plt
#import os # for reading/writing and directory operations

"""-----------------------------------------------------
Initialization: 
define lat/long bin grid, set read and write directories
--------------------------------------------------------"""

# extent of lat/long grid

# full manhattan (and a little brooklyn) extent
#bottom_left=[40.697089, -74.027397]
#top_left=[40.823067, -74.027397]
#bottom_right=[40.697089, -73.914240]
#top_right=[40.823067,-73.914240]

# wall street extent
#bottom_left=[40.701035, -74.022382]
#top_left=[40.715075, -74.022382]
#bottom_right=[40.701035, -74.000207]
#top_right=[40.715075, -74.000207]

# columbus circle extent
bottom_left=[40.765505, -73.985891]
top_left=[40.770369, -73.985891]
bottom_right=[40.765505,-73.978744]
top_right=[40.770369,-73.978744]

# set lat/long discretization
nx = 40 # number of bins in x (longitudinal) direction
ny = 50 # number of bins in y (latitudinal) direction

# the directory with the data
dir_base='./text_data/sub_sampling/'

# make the write directory 
write_dir='./data_columb_' + str(nx) + '_' + str(ny)
try: 
    os.makedirs(write_dir)
except OSError:
    if not os.path.isdir(write_dir):
       raise

# build cell node grid
x=np.linspace(bottom_left[1],bottom_right[1],nx) # defined on cell corners
y=np.linspace(bottom_left[0],top_left[0],ny) # defined on cell corners


"""----------------- count_taxi_file ----------------------- """
def count_taxi_file(f,nx,ny,x,y,there_is_a_header):
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
    TaxiCount=np.zeros((ny-1,nx-1)) # defined on cell centers

#   count number of lines, define the latitude, longitude arrays
    indx=0
    for line in f:
        indx=indx+1
    
    if there_is_a_header:
        indx = indx-1

    BigLat=np.zeros((indx,1))
    BigLon=np.zeros((indx,1))

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
           longi=float(line[5])
           lati=float(line[6])
           BigLat[indx,0]=lati
           BigLon[indx,0]=longi
           indx=indx+1
    
#   count the number of points within each lat/lon bin    
    for ix in range(nx-1):
        for iy in range(ny-1):
            ymin=y[iy] # lat limit 1
            ymax=y[iy+1] # lat limit 2
            xmin=x[ix] # long limit 1 
            xmax=x[ix+1] # long limit 2
            
            LocLat=BigLat[BigLon[:,0]>xmin]
            LocLon=BigLon[BigLon[:,0]>xmin]
    
            LocLat=LocLat[LocLon[:,0]<xmax]
            LocLon=LocLon[LocLon[:,0]<xmax]
    
            LocLon=LocLon[LocLat[:,0]>ymin]
            LocLat=LocLat[LocLat[:,0]>ymin]
    
            LocLon=LocLon[LocLat[:,0]<ymax]
            LocLat=LocLat[LocLat[:,0]<ymax]
    
            TaxiCount[iy,ix]=len(LocLon)
    return TaxiCount
"""-------------- END count_taxi_file ----------------------- """

"""----------------------------------------------------------
Processing: 
Loops over data files, counts number of taxi pickups within
each lat/lon bin
-------------------------------------------------------------"""
TaxiCount=np.zeros((ny-1,nx-1)) # zero array to add counts to
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
        TaxiCount2=count_taxi_file(fle,nx,ny,x,y,True) 

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
np.savetxt(write_dir+'/Nodes_x.txt', x, delimiter=',')
np.savetxt(write_dir+'/Nodes_y.txt', y, delimiter=',')

"""-------------------------- end of taxi_map_plot.py ----------------------"""
