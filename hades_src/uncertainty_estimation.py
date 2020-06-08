import numpy as num
import matplotlib.pyplot as plt

seq='ridgecrest'#'napa' or 'ridgecrest'

if seq=='napa':
   n=81
elif seq=='ridgecrest':
   n=320

data=num.load('locboot_real_'+seq+'.npy')
print(num.shape(data))
x_data_new=num.reshape(data[1:,0],(n,25))
y_data_new=num.reshape(data[1:,1],(n,25))
z_data_new=num.reshape(data[1:,2],(n,25))

x_data=[]
y_data=[]
z_data=[]
for i in range(n):
    x_data.append(num.std(x_data_new[i,:]))
    y_data.append(num.std(y_data_new[i,:]))
    z_data.append(num.std(z_data_new[i,:]))

x_data=num.sort(x_data)
y_data=num.sort(y_data)
z_data=num.sort(z_data)

c1='#4285F4'
c2='#EA4335'
c4='#34A853'
c3='#FBBC05'
n_bins=10

ax1=plt.subplot(1,1,1)
ax1.hist(x_data, n_bins, density=True, histtype='step', color=c1, cumulative=True)
#ax1.hist(refcat, n_bins, density=True, histtype='step', color=c2)
#ax1.set_xlim([0,30000])
#ax1.set_ylim([0,1])
#ax1.set_xticks(range(0,35000,5000))
#ax1.set_yticklabels(['%2.1f'%(i*0.1) for i in range(0,12,2)],fontsize='x-large')
#ax1.set_xticklabels(range(0,35,5),fontsize='x-large')
#ax1.set_xticks([])
#ax1.set_xlabel('Distance from HypoDD Locations (km)',fontsize='xx-large')
#ax1.set_ylabel('Fraction of events',fontsize='xx-large')
ax1.grid('on')

plt.savefig('histogram.eps')
plt.show()
