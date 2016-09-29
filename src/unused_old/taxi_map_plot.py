import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# the directory with the data to plot
dir_base='./data_man_80_100'
#dir_base='./data_110_120'
#dir_base='./data_wallst_110_120'
#dir_base='./data_columb_40_50'

# load the data 
TaxiCount=np.loadtxt(dir_base+'/TaxiCount.txt',delimiter=',')
MeanDist=np.loadtxt(dir_base+'/TaxiMeanDist.txt',delimiter=',')
x=np.loadtxt(dir_base+'/Nodes_x.txt',delimiter=',')
y=np.loadtxt(dir_base+'/Nodes_y.txt',delimiter=',')

# create grid info

# full meshgrid of nodes
[Xgrid,Ygrid]=np.meshgrid(x,y) 

# cell centers
xc = (x[0:x.size-1]+x[1:x.size])/2
yc = (y[0:y.size-1]+y[1:y.size])/2
[Xgridc,Ygridc]=np.meshgrid(xc,yc) 

# for color axis
minTax=2
maxTax=500#1000
minDist=0
maxDist=8


# smoothing functions
def two_D_average(window,RA):
    """  a two D average, loops over a 2D array and finds the mean of all points 
        within a radius of nodes. 

        input: 
            window    integer number of nearby nodes to average over
            RA        the twoD ndarray 
        output: 
            RAsmooth  the smoothed array

	nodes on edge of RA are not touched.     
    """
    n1 = RA.shape[0] # first dimension length
    n2 = RA.shape[1] # second dimension length

    NoTrouble = (n1-window) > window
    NoTrouble = NoTrouble * ((n2-window) > window)

    RAsmooth = RA

    if NoTrouble == 1: 
       for i1 in range(window,n1-window):
           for i2 in range(window,n2-window):
	        
		indx_1 = range(i1-window-1,i1+window+1)
		indx_2 = range(i2-window-1,i2+window+1)
		
		RAsmooth[i1,i2]=RA[indx_1,indx_2].mean()

    return(RAsmooth)


# pcolormesh plot
#Z = np.ma.masked_array(TaxiCount,mask = (TaxiCount<minTax))
#Z = np.ma.masked_array(Z,mask = (TaxiCount>maxTax))
Z = np.ma.masked_array(MeanDist,mask = (TaxiCount<minTax))
cm.hot.set_bad('black', alpha=None) 

#Z = TaxiCount * (TaxiCount > 3)#minTax) #*  (TaxiCount < maxTax)
#Z = MeanDist * (TaxiCount > 3)
#Z=two_D_average(1,Z)
#Z = np.log10(Z+.01)
#minTax=np.log10(minTax)
#maxTax=np.log10(maxTax)

pcol=plt.pcolormesh(Xgrid,Ygrid,Z,cmap=cm.hot,linewidth=0)
pcol.set_edgecolor('face')
#pcol=plt.pcolormesh(Xgrid,Ygrid,Z,cmap=cm.hot,linewidth=0)
#pcol.set_edgecolor('face')
#pcol=plt.pcolormesh(Xgrid,Ygrid,Z,cmap=cm.hot,linewidth=0)
#pcol.set_edgecolor('face')
#plt.clim(minTax,maxTax)
plt.clim(minDist,maxDist)
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
#plt.savefig('/home/chris/Desktop/foo.eps',figsize=(6,4)) 
