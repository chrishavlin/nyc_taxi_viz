import taxi_plotmod as tpm
import taxi_main as tm
import numpy as np
import matplotlib.pyplot as plt

# select directory and variable
data_dir='../data_products/three_month'
varname='dist_mi'

# output figure names
fig_name='dist_map'

# load binned data 
VarMean,VarCount,Varx,Vary=tm.read_gridded_file(data_dir,varname)   
VarMean=VarMean * (VarCount > 5)

# plot settings
log_plot=False
minval=0.1
maxval=10
x_dim=5
y_dim=7

ShowFig=False
SaveFig=False
figdpi=1200
savename=data_dir+'/'+fig_name

# plot it
tpm.plt_map(VarMean,minval,maxval,Varx,Vary,log_plot,ShowFig,SaveFig,savename,x_dim,y_dim,figdpi)


# calculate and plot mean trip distance vs latitude

Mean_vs_y = np.zeros((Vary.size,1))

for iLat in range(0,len(Vary)-1):
    x_var=VarMean[iLat,:]
    x_var=x_var[np.where(x_var>0.1)]
    Mean_vs_y[iLat] = np.mean(x_var)

fig=plt.figure()
plt.plot(Mean_vs_y,Vary)
plt.xlabel('Mean trip distance in miles')
plt.ylabel('Degree latitude')
axes = plt.gca()
axes.set_ylim([Vary.min(),Vary.max()])
axes.set_xlim([2,4.5])


plt.show()



    


