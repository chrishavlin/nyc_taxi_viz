"""---------------------------------------------------------------------------
Plotting routines

(c) Chris Havlin
Open source license?
----------------------------------------------------------------------------"""

"""--------------
Import libraries:
-----------------"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

"""---------
Classes
------------"""
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
        self.mean=np.zeros((self.time_be.size-1,1))
        self.med=np.zeros((self.time_be.size-1,1))
        self.std=np.zeros((self.time_be.size-1,1))
        self.N=np.zeros((self.time_be.size-1,1))

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

"""---------
Functions
------------"""
def map_proc(VarBig,Var_list,bin_varname,bin_min,bin_max,Pickup,nx,ny,
            bottom_left=[40.697089, -74.027397],top_left=[40.823067, -74.027397],
            bottom_right=[40.697089, -73.914240],top_right=[40.823067,-73.914240]):
   ##############################################################################
   # Input Variables:
   # VarBig
   #    The 2-D array of all variables
   # Var_list
   #    The list identifying the columns of VarBig
   # bin_varname
   #    The name of the variable to bin
   # bin_min, bin_max
   #    The min/max value of the variable to bin   
   # Pickup
   #    Boolean true/false. If true, the variable of interest (bin_varname)
   #    will be binned by the pick up location. If false, the drop off location. 
   # nx,ny
   #    Number of bin-nodes in x and y (i.e., number of bins in x is nx-1)
   # bottom_left, top_left, etc.    
   #    Grid dimensions: lat/lon points defining the corners of a rectangle
   #
   # Some useful lat/lon points	    
   # Manhattan:
   #        bottom_left=[40.697089, -74.027397],top_left=[40.823067, -74.027397],
   #        bottom_right=[40.697089, -73.914240],top_right=[40.823067,-73.914240]
   # Wall Street:
   #        bottom_left=[40.701035, -74.022382],top_left=[40.715075, -74.022382],
   #        bottom_right=[40.701035, -74.000207],top_right=[40.715075,-74.000207]
   # Columbus Circle:
   #        bottom_left=[40.765505, -73.985891],top_left=[40.770369, -73.985891],
   #        bottom_right=[40.765505, -73.978744],top_right=[40.770369,-73.978744]
   ##############################################################################

    # build cell node grid
    x=np.linspace(bottom_left[1],bottom_right[1],nx) # defined on cell corners
    y=np.linspace(bottom_left[0],top_left[0],ny) # defined on cell corners

    # cell center grid
    xc=(x[0:x.size-1]+x[1:x.size])/2
    yc=(y[0:y.size-1]+y[1:y.size])/2

    # create an empty TaxiCount array 
    TotCount=np.zeros((ny-1,nx-1)) # defined on cell centers
    MeanClr=np.zeros((ny-1,nx-1)) # defined on cell centers

    # pull out lat/lon and variable lists depending on pickup or dropoff
    if Pickup:
       BigLat=VarBig[:,Var_list.index('pickup_lat')] 
       BigLon=VarBig[:,Var_list.index('pickup_lon')] 
    else:
       BigLat=VarBig[:,Var_list.index('drop_lat')] 
       BigLon=VarBig[:,Var_list.index('drop_lon')] 

    ColorVar= VarBig[:,Var_list.index(bin_varname)]
    indc=1

#    # count the number of points within each lat/lon bin
#    for ix in range(nx-1):
#        for iy in range(ny-1):
#            ymin=y[iy] # lat limit 1
#            ymax=y[iy+1] # lat limit 2
#            xmin=x[ix] # long limit 1 
#            xmax=x[ix+1] # long limit 2
#            
#            LocClr=ColorVar[np.logical_and(BigLon>xmin,BigLon<xmax)]
#            LocLat=BigLat[np.logical_and(BigLon>xmin,BigLon<xmax)]
#            
#            LocClr=LocClr[np.logical_and(LocLat>ymin,LocLat<ymax)]
#            LocClr=LocClr[np.logical_and(LocClr>bin_min,LocClr<bin_max)]
#
#            if len(LocClr)>0:
#               TotCount[iy,ix]=len(LocClr)
#               MeanClr[iy,ix]=LocClr.mean()
#               #print indc,'of',(nx*ny),'val:',TotCount[iy,ix]
#            #else:
#            #   #print indc,'of',(nx*ny),'val: empty'
#            indc=indc+1

    # count, but loop over variable
    nV = ColorVar.size
    dx = abs(xc[1]-xc[0])
    dy = abs(yc[1]-yc[0])
    prevprog=0
    print 'Binning',bin_varname,', with Nvar=',nV
    for iV in range(nV-1):
        prog= round(float(iV) / float(nV-1) * 100)
        
        if prog % 5 == 0 and prog != prevprog:
           print '    ',int(prog),'% complete ...'
           prevprog=prog

        if ColorVar[iV]>bin_min and ColorVar[iV]<bin_max:
           xi=BigLon[iV] 
           yi=BigLat[iV]
           i_y=np.where(abs(yi-yc)<dy/2.0)
           i_x=np.where(abs(xi-xc)<dx/2.0)
           
           if i_y[0].size==1 and i_x[0].size==1:
              TotCount[i_y[0],i_x[0]]=TotCount[i_y[0],i_x[0]]+1
              MeanClr[i_y[0],i_x[0]]=MeanClr[i_y[0],i_x[0]]+ColorVar[iV]

    non0=np.where(TotCount>0)
    if non0[0].size>0:
       MeanClr[non0]=np.divide(MeanClr[non0],TotCount[non0])

    print 'Completed binning of ',bin_varname

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

def plt_map(Zvar,minZ,maxZ,x,y,LogPlot=False,ShowFig=True,SaveFig=False,savename=' ',
            dim_x_in=4,dim_y_in=6,figdpi=180):

    fig=plt.figure()
    fig.set_size_inches(dim_x_in,dim_y_in,forward=True)

    if LogPlot==True:
       Zvar = np.log10(Zvar)
       minZ = np.log10(minZ)
       maxZ = np.log10(maxZ)

    # mask the variable
    Z = Zvar
    
    # cell edges
    [Xgrid,Ygrid]=np.meshgrid(x,y) 

    # get cell centers and make the meshgrid
    xc = (x[0:x.size-1]+x[1:x.size])/2
    yc = (y[0:y.size-1]+y[1:y.size])/2
    [Xgridc,Ygridc]=np.meshgrid(xc,yc)

    # now plot it!
    pcol=plt.pcolormesh(Xgrid,Ygrid,Z,cmap=cm.hot,linewidth=0)
    pcol.set_edgecolor('face')
    plt.clim(minZ,maxZ)
    plt.colorbar()
    
    if ShowFig:
      print 'close figure to continue...'
      plt.show()

    if SaveFig:
      fig.savefig(savename+'.png',format='png', dpi=figdpi)
