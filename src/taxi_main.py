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
    """ class for variables binned by time """
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
         
""" END OF FUNCTIONS AND CLASSES """
""" PLOTTING FUNCTIONS BELOW """
def map_proc(VarBig,Var_list,bin_varname,bin_min,bin_max,Pickup,nx,ny,
            bottom_left=[40.697089, -74.027397],top_left=[40.823067, -74.027397],
            bottom_right=[40.697089, -73.914240],top_right=[40.823067,-73.914240]):
    # default is full manhattan

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

    # build cell node grid
    x=np.linspace(bottom_left[1],bottom_right[1],nx) # defined on cell corners
    y=np.linspace(bottom_left[0],top_left[0],ny) # defined on cell corners

    # create an empty TaxiCount array 
    TotCount=np.zeros((ny-1,nx-1)) # defined on cell centers
    MeanClr=np.zeros((ny-1,nx-1)) # defined on cell centers

    # pull out lat/lon and variable lists depending on pickup or dropoff
    if Pickup:
       BigLat=VarBig[:,8] 
       BigLon=VarBig[:,7] 
    else:
       BigLat=VarBig[:,10] 
       BigLon=VarBig[:,9] 

    ColorVar= VarBig[:,Var_list.index(bin_varname)]
    indc=1

    # count the number of points within each lat/lon bin
    for ix in range(nx-1):
        for iy in range(ny-1):
            ymin=y[iy] # lat limit 1
            ymax=y[iy+1] # lat limit 2
            xmin=x[ix] # long limit 1 
            xmax=x[ix+1] # long limit 2
            
            LocClr=ColorVar[np.logical_and(BigLon>xmin,BigLon<xmax)]
            LocLat=BigLat[np.logical_and(BigLon>xmin,BigLon<xmax)]
            
            LocClr=LocClr[np.logical_and(LocLat>ymin,LocLat<ymax)]
            LocClr=LocClr[np.logical_and(LocClr>bin_min,LocClr<bin_max)]

#            LocClr=ColorVar[BigLon[:]>xmin]
#            LocLat=BigLat[BigLon[:]>xmin]
#            LocLon=BigLon[BigLon[:]>xmin]
#
#            LocClr=LocClr[LocLon[:]<xmin]
#            LocLat=LocLat[LocLon[:]<xmax]
#            LocLon=LocLon[LocLon[:]<xmax]
#
#            LocClr=LocClr[LocLat[:]>ymin]
#            LocLat=LocLat[LocLat[:]>ymin]
#
#            LocClr=LocClr[LocLat[:]<ymax]

            #LocClr=LocClr[LocClR[:]>bin_min]
            #LocClr=LocClr[LocClR[:]<bin_max]

            if len(LocClr)>0:
               TotCount[iy,ix]=len(LocClr)
               MeanClr[iy,ix]=LocClr.mean()
               print indc,'of',(nx*ny),'val:',TotCount[iy,ix]
            else:
               print indc,'of',(nx*ny),'val: empty'
            indc=indc+1

    return TotCount,MeanClr,x,y


def plt_two_d_histogram(bin_varname,VarMin,VarMax,time_b,VarBig,Var_list):
    
    bin_inst=binned_variable(bin_varname,time_b,VarMin,VarMax)
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

def plt_map(Zvar,minZ,maxZ,x,y):
    plt.figure()

    # mask the variable
    Z = Zvar
    #Z = np.ma.masked_array(Zvar,mask = (Zvar>minZ))
    #Z = np.ma.masked_array(Z,mask = (Z<maxZ))
    #cm.hot.set_bad('black', alpha=None)
    
    # cell edges
    [Xgrid,Ygrid]=np.meshgrid(x,y) 

    # cell centers
    xc = (x[0:x.size-1]+x[1:x.size])/2
    yc = (y[0:y.size-1]+y[1:y.size])/2
    [Xgridc,Ygridc]=np.meshgrid(xc,yc)

    pcol=plt.pcolormesh(Xgrid,Ygrid,Z,cmap=cm.hot,linewidth=0)
    pcol.set_edgecolor('face')
    plt.clim(minZ,maxZ)
    plt.colorbar()
    
    # scatter plot
    #nconts=300
    #Z=np.multiply(TaxiCount,TaxiCount>minTax)
    #Z=np.multiply(Z,Z<maxTax)
    #plt.scatter(Xgridc,Ygridc,c=Z,s=100,cmap=cm.hot,edgecolors=None,linewidths=0)
    #plt.clim(minTax,maxTax)
    #plt.colorbar()
    
    # contour plot
    #plt.subplot(1,2,1)
    #nconts=400
    #plt.contourf(Xgridc,Ygridc,TaxiCount,nconts,cmap=cm.hot)
    #plt.contourf(Xgridc,Ygridc,TaxiCount,nconts,cmap=cm.hot)
    #plt.clim(minTax,maxTax)
    #plt.colorbar()
    
    # log plot
    #nconts=300
    #plt.subplot(1,2,2)
    #plt.contourf(Xgridc,Ygridc,np.log10(TaxiCount),nconts,cmap=cm.hot)
    #plt.clim(np.log10(minTax),np.log10(maxTax))
    #plt.colorbar()
    
    plt.show()

""" END OF PLOTTING FUNCTIONS """

if __name__ == '__main__':
    # the directory with the data
    #dir_base='../data_full_textdata/'
    dir_base='../data_full_textdata/sub_sampling_16/'

    # read in all the data! 
    VarBig,Var_list=read_taxi_files(dir_base)
    
    """--------------------------------------------------"""
    #Var_list=['pickup_time_hr','dist_mi','speed_mph','psgger','fare',
    #          'tips','payment_type','pickup_lon','pickup_lat','drop_lon',
    #          'drop_lat']
     

    FareCount,FareMean,Farex,Farey=map_proc(VarBig,Var_list,'fare',1,100,'True',50,80)

    plt_map(FareCount,1,FareCount.max(),Farex,Farey)

    bin_varname = 'speed_mph' # the variable to bin 
    time_b = np.linspace(0,24,150) # time bin edges
    plt_two_d_histogram(bin_varname,1,60,time_b,VarBig,Var_list)
  

#    bin_varname = 'psgger' # the variable to bin 

