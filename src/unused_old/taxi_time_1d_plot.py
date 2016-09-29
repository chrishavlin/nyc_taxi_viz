import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# the directory with the data to plot
dir_base='./data_timeN__96'

# load the data 
TaxiCount=np.loadtxt(dir_base+'/TaxiCount.txt',delimiter=',')
TaxiCount=TaxiCount/31.0
Time=np.loadtxt(dir_base+'/Time.txt',delimiter=',')

# create grid info

# cell centers
Timec = (Time[0:Time.size-1]+Time[1:Time.size])/2
dt=Time[1]-Time[0]

# for color axis
plt.subplot(2,1,2)
plt.plot(Timec,TaxiCount/(dt*60))
plt.xlim(0,24)
plt.ylabel('Pick ups per minute')
plt.xlabel('Time of day (24-hour time)')

plt.subplot(2,1,1)
plt.plot(Timec,TaxiCount)
plt.xlim(0,24)
plt.ylabel('Pick ups in 15 minute bins (daily average)')
plt.xlabel('Time of day (24-hour time)')

plt.show()
#plt.savefig('/home/chris/Desktop/foo.eps',figsize=(6,4)) 
