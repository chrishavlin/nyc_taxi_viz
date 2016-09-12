"""---------------------------------------------------------------------------
Processes the csv taxi files. 

(c) Chris Havlin
Open source license?
----------------------------------------------------------------------------"""

"""--------------
Import libraries:
-----------------"""
import numpy as np
import time,os # for 
import matplotlib.pyplot as plt
from matplotlib import cm

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
    Vars=np.zeros((indx,7))

#   Go back to start of file, loop again to read variables
    f.seek(0)

    if there_is_a_header:
        headerline=f.readline()

    indx=0
    Var_list=['pickup_time_hr','dist_mi','speed_mph','psgger','fare',
              'tips','payment_type']

#   loop over lines, store variables
    for line in f:
        line = line.rstrip()
        line = line.split(',')
        if len(line) == 19:
	 # pickup time
	   pickup=line[1] # datetime string
	   pickup=pickup.split() # remove the space
	   pickup=pickup[1] # take time only (date is pickup[0])
	   t_hms=pickup.split(':')
	   pickup_time=float(t_hms[0])+float(t_hms[1])/60+float(t_hms[2])/3600

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

         # store in Vars 
           Vars[indx,0]=pickup_time
           Vars[indx,1]=dist
           Vars[indx,2]=speed
           Vars[indx,3]=psgger
           Vars[indx,4]=fare
           Vars[indx,5]=tips
           Vars[indx,6]=pay_type

           indx=indx+1

    return Vars,Var_list
"""-------------- END count_taxi_file ----------------------- """

def read_taxi_files(dir_base):

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

class binned_variable(object):
    def __init__(self,name,time_b,varmin=0,varmax=1):
        self.varname = name # variable name
        self.time_be = time_b # time bins, edge values
        self.time_bc = (time_b[0:time_b.size-1]+time_b[1:time_b.size])/2.0 # time bin center values
        self.varmin=varmin
        self.varmax=varmax
    
    def bin_the_values(self,VarBig,Var_list):
        # initialize the output 
        self.mean=np.zeros((time_b.size-1,1))
        self.med=np.zeros((time_b.size-1,1))
        self.std=np.zeros((time_b.size-1,1))
        self.N=np.zeros((time_b.size-1,1))

        # extract variable of interest from full dataset
        var0 = VarBig[:,Var_list.index(self.varname)]
        time0 = VarBig[:,0]

        # calculate standard deviation,mean,median and number of obs in each bin
        for i_bin in range(self.time_bc.size):
            var2bin = var0[time0>self.time_be[i_bin]]
            time = time0[time0>self.time_be[i_bin]]
        
            var2bin = var2bin[time<self.time_be[i_bin+1]]
        
            var2bin = var2bin[var2bin>self.varmin]
            var2bin = var2bin[var2bin<self.varmax]
    
            self.mean[i_bin]=var2bin.mean()
            self.med[i_bin]=np.median(var2bin)
            self.std[i_bin]=var2bin.std()
	    self.N[i_bin]=var2bin.size

    def bin_hist(self,VarBig,Var_list,bin_edge1,bin_edge2):
        # extract variable of interest from full dataset
        var0 = VarBig[:,Var_list.index(self.varname)]
        time0 = VarBig[:,0]

	bin_var=var0[time0>bin_edge1]
	time=time0[time0>bin_edge1]

	bin_var=bin_var[time<bin_edge2]
        bin_var = bin_var[bin_var>self.varmin]
        bin_var = bin_var[bin_var<self.varmax]

	self.hist_id=bin_var
	self.hist_id_N=bin_var.size
	self.hist_id_mean=bin_var.mean()
	self.hist_id_med=np.median(bin_var)
         

    
if __name__ == '__main__':
    # the directory with the data
    #dir_base='./text_data/sub_sampling_16/'
    #dir_base='./text_data/single_file/'
    dir_base='./text_data/'

    VarBig,Var_list=read_taxi_files(dir_base)

    
    """--------------------------------------------------"""
    # Var_list=['pickup_time_hr','dist_mi','speed_mph','psgger','fare',
    #           'tips','payment_type']
     

    bin_varname = 'speed_mph' # the variable to bin 
#    bin_varname = 'psgger' # the variable to bin 
    time_b = np.linspace(0,24,150) # time bin edges
    
    bin_inst=binned_variable(bin_varname,time_b,1,60)
    bin_inst.bin_the_values(VarBig,Var_list)

    # all values, limited by min,max

    max_unbin=bin_inst.varmax
    min_unbin=bin_inst.varmin

    unbinned_value = VarBig[:,Var_list.index(bin_varname)]
    time = VarBig[:,0]
    
    time = time[unbinned_value>=min_unbin]
    unbinned_value = unbinned_value[unbinned_value>=min_unbin]
    
    time = time[unbinned_value<=max_unbin]
    unbinned_value = unbinned_value[unbinned_value<=max_unbin]


    clr = (0.3,0.3,0.3)
    
    plt.subplot(1,3,1)
    plt.hist2d(time,unbinned_value,(48,30),cmap=cm.hot)
    plt.colorbar()
    plt.plot(bin_inst.time_bc,bin_inst.mean,color=(0.6,0.6,0.6),linewidth=3)
    plt.plot(bin_inst.time_bc,bin_inst.med,color=clr,linewidth=3)
    plt.xlim([0,24])
    plt.ylim([min_unbin,max_unbin])
    plt.ylim([min_unbin,50])
    plt.ylabel(bin_varname)
    plt.xlabel('time of day [24-hr]')
    
    plt.subplot(1,3,2)
    bin_inst.bin_hist(VarBig,Var_list,5,6)
    LAB='5-6,'+str(bin_inst.hist_id_N) + ',' + str(round(bin_inst.hist_id_med,1))
    LAB = LAB + ',' + str(round(bin_inst.hist_id_mean,1))
    plt.hist(bin_inst.hist_id,bins=50,histtype='step',normed=True,label=LAB)

    bin_inst.bin_hist(VarBig,Var_list,9,18)
    LAB='9-18,' +str(bin_inst.hist_id_N) + ',' + str(round(bin_inst.hist_id_med,1))
    LAB = LAB + ',' + str(round(bin_inst.hist_id_mean,1))
    plt.hist(bin_inst.hist_id,bins=50,histtype='step',normed=True,label=LAB)

    bin_inst.bin_hist(VarBig,Var_list,20,24)
    LAB='20-24,'+str(bin_inst.hist_id_N) + ',' + str(round(bin_inst.hist_id_med,1))
    LAB = LAB + ',' + str(round(bin_inst.hist_id_mean,1))
    plt.hist(bin_inst.hist_id,bins=50,histtype='step',normed=True,label=LAB)
    plt.ylabel('bin N / tot N')

    plt.xlabel(bin_varname)
    plt.ylim([0,0.15])
    plt.legend()
    
    plt.subplot(1,3,3)
    plt.plot(bin_inst.time_bc,2*bin_inst.std)
    plt.xlabel('time of day [24-hr]')
    plt.ylabel('2-standard deviations of ' + bin_varname)
    plt.xlim([0,24])

    #LAB=str(unbinned_value.size) + ',' + str(round(np.median(unbinned_value),1))
    
    plt.figure()
    plt.plot(bin_inst.N,bin_inst.med)
    plt.scatter(bin_inst.N,bin_inst.med,c=bin_inst.time_bc)
    plt.colorbar()
    plt.xlabel('N')
    plt.ylabel('speed [mph]')
    plt.show()
    
    
    #np.savetxt(write_dir+'/Vars.txt', Vars, delimiter=',')
