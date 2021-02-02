import matplotlib.pyplot as plt
import numpy as num
km=1000

def read_catalogue(cat_file, pos=-1):
    with open(cat_file, 'r') as f:
        catlocs=[]
        for line in f:
            toks=line.split()
            depth=eval(toks[pos])*km
            catlocs.append(depth)
    catloc=num.array(catlocs)
    catloc=num.sort(catloc)
    return catloc

refcat=read_catalogue('./napa_data/napa_refcat_full_DD.dat')
#cat_WBS=read_catalogue('outcat_ridgecrest_WBS.txt')
#cat_WMF=read_catalogue('outcat_ridgecrest_WMF.txt')
cat=read_catalogue('outcat_napa.txt')

c1='#4285F4'
c2='#EA4335'
c4='#34A853'
c3='#FBBC05'
n_bins=10


ax1=plt.subplot(1,1,1)
ax1.hist(cat, n_bins, density=True, histtype='step', color=c1)
ax1.hist(refcat, n_bins, density=True, histtype='step', color=c2)
ax1.set_xlim([0,30000])
#ax1.set_ylim([0,1])
ax1.set_xticks(range(0,35000,5000))
#ax1.set_yticklabels(['%2.1f'%(i*0.1) for i in range(0,12,2)],fontsize='x-large')
#ax1.set_xticklabels(range(0,35,5),fontsize='x-large')
ax1.set_xticks([])
#ax1.set_xlabel('Distance from HypoDD Locations (km)',fontsize='xx-large')
#ax1.set_ylabel('Fraction of events',fontsize='xx-large')
ax1.grid('on')

plt.savefig('histogram.eps')
plt.show()
