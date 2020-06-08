import molquake
import matplotlib.pyplot as plt
import numpy as num


output_path='/Users/francesco/Desktop/Single_Station_Location/MOLQUAKE_V_0.3/real_data'
data_path='/Users/francesco/Desktop/Single_Station_Location/MOLQUAKE_V_0.3/real_data'
datafiles=['napa_refcat_4.dat','napa_picks.dat']
out_catalogue='outcat.txt'
datafile='vp6800'
Vp=6800 #velocities in m/s !!!!!
Vps=1.73
Vs=Vp/Vps
nref=4

dataobj=molquake.molquake_dat(data_path,datafiles)
locobj=molquake.molquake_loc(Vp,Vs)
references=dataobj.references
nref=num.size(references[:,0])
locobj.location(dataobj.data, references)
locobj.catalogue_creation(dataobj.evids, dataobj.origin, nref, out_catalogue, output_path)
#multiloc=locobj.uncertainty_estimation(dataobj.data, references, 25, [6500,7500], 1.73)
#print(num.shape(multiloc))
#locboot=num.array([0.,0.,0])
#for i in range(20):
#      locboot=num.vstack((locboot,multiloc[i,:,:]))
#num.save('locboot_real.npy',locboot)
#num.save('locreal.npy',dataobj.references)


nref,mref=num.shape(references)
err_svd=num.sort(num.sqrt((dataobj.references[:,0]-locobj.locations[:,0])**2+(dataobj.references[:,1]-locobj.locations[:,1])**2+(dataobj.references[:,2]-locobj.locations[:,2])**2))

n_bins=30
#c1='#FF6347'
#c2='#4682B4'
#c3='#99e699'
msize=30
c1='#4285F4'
c2='#EA4335'
c4='#34A853'
c3='#FBBC05'

fig=plt.figure(figsize=(10,10))
ax1=plt.subplot(2,2,1)
ax1.scatter(locobj.locations[:,0],locobj.locations[:,1], s=msize, c=c1)
ax1.scatter(dataobj.references[:,0],dataobj.references[:,1],s=msize, c=c2)
ax1.scatter(references[:,0],references[:,1], s=msize, c=c3)
ax1.set_xlim([-10000,10000])
ax1.set_ylim([-10000,10000])
ax1.set_xticks(range(-10000,12000,2000))
ax1.set_yticks(range(-10000,12000,2000))
ax1.set_xticklabels(range(-10,12,2))
ax1.set_yticklabels(range(-10,12,2))
ax1.set_xlabel('X(km)')
ax1.set_ylabel('Y(km)')
ax1.grid('on')
ax1.set_aspect('equal')

ax2=plt.subplot(2,2,2)
ax2.scatter(locobj.locations[:,2],locobj.locations[:,1],s=msize, c=c1)
ax2.scatter(dataobj.references[:,2],dataobj.references[:,1],s=msize, c=c2)
ax2.scatter(references[:,2],references[:,1],s=msize, c=c3)
ax2.set_xlim([0,20000])
ax2.set_ylim([-10000,10000])
ax2.set_xticks(range(0,22000,2000))
ax2.set_yticks(range(-10000,12000,2000))
ax2.set_xticklabels(range(0,22,2))
ax2.set_yticklabels(range(-10,12,2))
ax2.set_xlabel('Z(km)')
ax2.set_ylabel('Y(km)')
ax2.set_aspect('equal')
ax2.grid('on')

ax3=plt.subplot(2,2,3)
ax3.scatter(locobj.locations[:,0],locobj.locations[:,2],s=msize, c=c1)
ax3.scatter(dataobj.references[:,0],dataobj.references[:,2],s=msize, c=c2)
ax3.scatter(references[:,0],references[:,2],s=msize, c=c3)
ax3.set_ylim([0,20000])
ax3.set_xlim([-10000,10000])
ax3.set_yticks(range(0,22000,2000))
ax3.set_xticks(range(-10000,12000,2000))
ax3.set_yticklabels(range(0,22,2))
ax3.set_xticklabels(range(-10,12,2))
ax3.set_xlabel('X(km)')
ax3.set_ylabel('Z(km)')
ax3.invert_yaxis()
ax3.set_aspect('equal')
ax3.grid('on')

ax4=plt.subplot(2,2,4)
#ax4.hist(err_svd, n_bins, density=True, histtype='step', cumulative=True, color=c4)
ax4.set_xlim([0,10000])
ax4.set_ylim([0,1])
ax4.set_xlabel('Distance from True Location (m)')
ax4.set_ylabel('Fraction of events')
ax4.grid('on')

plt.savefig(datafile+'_'+str(nref)+'.eps')
plt.show()

plt.plot(multiloc[:,:,0],multiloc[:,:,1],'ro')
plt.plot(dataobj.references[:,0],dataobj.references[:,1],'b+')
plt.xlim([-10000,10000])
plt.ylim([-15000,15000])
plt.grid('on')
plt.show()

plt.plot(multiloc[:,:,0],multiloc[:,:,2],'ro')
plt.plot(dataobj.references[:,0],dataobj.references[:,2],'b+')
plt.xlim([-15000,15000])
plt.ylim([-15000,15000])
plt.grid('on')
plt.show()

plt.plot(multiloc[:,:,1],multiloc[:,:,2],'ro')
plt.plot(dataobj.references[:,1],dataobj.references[:,2],'b+')
plt.xlim([-15000,15000])
plt.ylim([-15000,15000])
plt.grid('on')
plt.show()
