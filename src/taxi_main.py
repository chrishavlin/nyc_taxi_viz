"""---------------------------------------------------------------------------
Processes the csv taxi files. 

(c) Chris Havlin
Open source license?
----------------------------------------------------------------------------"""

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
def read_all_variables(f,there_is_a_header):
    """ 
     reads in the data
    
     input: 
       f                   file object
       there_is_a_header   logical flag 
       
     output: 
       Vars
    """

#   count number of lines
    indx=0
    for line in f:
        indx=indx+1
    if there_is_a_header:
        indx = indx-1

#   Initizialize Variable Array
    Vars=np.zeros((indx,11))

#   Go back to start of file, loop again to read variables
    f.seek(0)

    if there_is_a_header:
        headerline=f.readline()

    indx=0
    Var_list=['pickup_time_hr','dist_mi','speed_mph','psgger','fare',
              'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
              'drop_lat']

#   loop over lines, store variables
    for line in f:
        line = line.rstrip()
        line = line.split(',')
        if len(line) == 19:
	 # pickup time
	   pickup=line[1].split()[1] # splits date-time string at the space, takes
	                             # the element that is the time string
	   t_hms=pickup.split(':') # splits into hours, min, sec
	   pickup_time=float(t_hms[0])+float(t_hms[1])/60+float(t_hms[2])/3600

	 # pickup date
	   date=line[1].split()[0] # the date string, "yyyy-mm-dd"

         # dropoff time
	   drop=line[2] # datetime string
	   drop=drop.split() # remove the space
	   drop=drop[1] # take time only (date is drop[0])
	   t_hms=drop.split(':')
	   drop_time=float(t_hms[0])+float(t_hms[1])/60+float(t_hms[2])/3600

         # distance travelled [mi]
           dist=float(line[4]) # distance travelled [mi]

	 # average speed
	   if (drop_time-pickup_time) > 0:
               speed=dist / (drop_time - pickup_time) # [mi/hr]
	   else:
	       speed=0
           
	 # other variables
	   psgger=float(line[3]) # passenger count
	   fare=float(line[12]) # base fare ammount [$]
	   tips=float(line[15]) # tip amount [$]
	   pay_type=float(line[11]) # payment type
	   pickup_lat=float(line[6])
           pickup_lon=float(line[5])
	   drop_lat=float(line[10])
           drop_lon=float(line[9])

         # store in Vars 
           Vars[indx,0]=pickup_time
           Vars[indx,1]=dist
           Vars[indx,2]=speed
           Vars[indx,3]=psgger
           Vars[indx,4]=fare
           Vars[indx,5]=tips
           Vars[indx,6]=pay_type
           Vars[indx,7]=pickup_lon
           Vars[indx,8]=pickup_lat
           Vars[indx,9]=drop_lon
           Vars[indx,10]=drop_lat

           indx=indx+1

    return Vars,Var_list

def read_taxi_files(dir_base):
    """ reads in individual taxi file, stores variables """

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
            VarChunk,Var_list=read_all_variables(fle,True) 
    
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

#def write_out_processed_file():   

#def read_in_processed_file():

""" END OF FUNCTIONS """

if __name__ == '__main__':
    # the directory with the data
    dir_base='../data_full_textdata/'
    #dir_base='../data_full_textdata/sub_sampling_16/'
    #dir_base='../full_csv_files/'

    # read in all the data! 
    VarBig,Var_list=read_taxi_files(dir_base)
    
    """--------------------------------------------------"""
    #Var_list=['pickup_time_hr','dist_mi','speed_mph','psgger','fare',
    #          'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
    #          'drop_lat']

    DistCount,DistMean,Distx,Disty=tpm.map_proc(VarBig,Var_list,'dist_mi',0.1,30,'True',600,700)

    DistMean = DistMean * (DistCount > 10)
    tpm.plt_map(DistMean,1,10,Distx,Disty,False)
    tpm.plt_map(DistCount,1,1000,Distx,Disty,True)

    bin_varname = 'speed_mph' # the variable to bin 
    time_b = np.linspace(0,24,10) # time bin edges
    tpm.plt_two_d_histogram(bin_varname,1,60,time_b,VarBig,Var_list)
  
