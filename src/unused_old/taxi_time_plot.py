import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import axis
from pylab import rcParams
rcParams['figure.figsize'] = 7.5, 5.5

# the directory with the data to plot
dir_base='./data_time_nt_48_daily'

# load the data 
TaxiCount=np.loadtxt(dir_base+'/TaxiCount.txt',delimiter=',')
Time_Bins=np.loadtxt(dir_base+'/Time.txt',delimiter=',')
Days=np.loadtxt(dir_base+'/Days.txt',delimiter=',')
Day_Bins=np.insert(Days,0,0)

# shift values so that midnight to 5am appears on previous day
TaxiShift=np.zeros(TaxiCount.shape)
Time_Shift=Time_Bins + 5
Time_c = (Time_Bins[0:Time_Bins.size-1] + Time_Bins[1:Time_Bins.size])/2
Time_c_Sh=Time_c+5

for iday0 in range(0,len(Days)):
    iday1 = iday0+1
    iday1 = iday1 * (iday1 <= len(Days)-1) + (0) * (iday1> len(Days)-1)

    night_rides=TaxiCount[Time_c<5,iday1]
    day_rides = TaxiCount[Time_c>=5,iday0]

    if iday1==0:
       night_rides[:]=-1

    TaxiShift[:,iday0]=np.concatenate((day_rides,night_rides))

Z = np.ma.masked_array(TaxiShift,mask = (TaxiShift<0))

# adjust tick locations and labels
fullrange=range(5,30,1)
part2=range(1,6,1)
part1=range(5,25,1)

# calculate some mean values
MeanTaxi = np.mean(TaxiShift,1) # total mean
DayMean = np.zeros((7,MeanTaxi.size))
DayCount = np.zeros((7,1))
DayOfWeek=0
for iday in range(0,len(Days)): 
    DayCount[DayOfWeek]=DayCount[DayOfWeek]+1
    DayMean[DayOfWeek,:]=DayMean[DayOfWeek,:]+TaxiShift[:,iday]

    DayOfWeek = DayOfWeek+1
    if DayOfWeek > 6:
       DayOfWeek = 6
    
for iday in range(0,7):
    DayMean[iday,:]=DayMean[iday,:] / DayCount[iday]
    DayMean[iday,:]=DayMean[iday,:]/DayMean[iday,:].max()

print DayMean.shape,DayCount.shape    

# set masked value color
cm.rainbow.set_bad('black', alpha=None)

# plot
plt.figure()
#plt.subplot(1,2,1)
#plt.pcolormesh(Time_Shift,Day_Bins,Z.transpose(),cmap=cm.rainbow)
##plt.pcolormesh(Time_Shift,Day_Bins,Z.transpose(),cmap=cm.rainbow)
##plt.pcolormesh(Time_Shift,Day_Bins,Z.transpose(),cmap=cm.rainbow)
#plt.xticks(range(5,30,1),part1+part2)
#plt.yticks(np.arange(0.5,33.5,1),range(1,33,1))
#plt.xlim(5,24+5)
#plt.ylim(0,Day_Bins.max())
#plt.colorbar()
#
##plt.savefig('/home/chris/Desktop/foo.eps')
##plt.show()



#TaxiCount


WeekDay=['Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday']
Clrs=[(0,0,0),(0,0,1),(0,0,1),(0,0,1),(0,1,0),(1,0,0),(0,1,0)]
#plt.subplot(1,2,2)
#plt.plot(Time_c,MeanTaxi/MeanTaxi.max())
for iday in range(0,7):
    scl=float(iday)/6
    lab=WeekDay[iday]
    #plt.plot(Time_c_Sh,DayMean[iday,:],color=[(scl),0,(1-scl)],label=lab)
    clr = Clrs[iday][:]
    plt.plot(Time_c_Sh,DayMean[iday,:],label=lab,color=clr)


plt.legend(loc=8)
plt.xlim(5,24+5)
plt.xticks(range(5,30,1),part1+part2)
plt.show()
