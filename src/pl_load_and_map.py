import taxi_plotmod as tpm
import taxi_main as tm
import numpy as np

# select directory and variable
data_dir='../data_products/three_month'
varname='speed_mph'

# output figure names
fig_name='speed_map'

# load binned data 
VarMean,VarCount,Varx,Vary=tm.read_gridded_file(data_dir,varname)   

# plot settings
log_plot=False
minval=1
maxval=60#VarCount.max()
x_dim=5
y_dim=7


ShowFig=False
SaveFig=True
savename=data_dir+'/'+fig_name
figdpi=580
# def plt_map(Zvar,minZ,maxZ,x,y,LogPlot=False,ShowFig=True,SaveFig=False,savename=' ',
# 257             dim_x_in=4,dim_y_in=6,figdpi=180):

# plot it
tpm.plt_map(VarMean,minval,maxval,Varx,Vary,log_plot,ShowFig,SaveFig,savename,x_dim,y_dim,figdpi)

