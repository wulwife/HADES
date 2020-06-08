import molquake
import matplotlib.pyplot as plt
import numpy as num

output_path='/Users/francesco/Desktop/Single_Station_Location/MOLQUAKE_V_0.3/napa_data'
data_path='/Users/francesco/Desktop/Single_Station_Location/MOLQUAKE_V_0.3/napa_data'
datafiles=['napa_refcat_16_DD.dat','napa_picks_16.dat']
catalogue_file='napa_refcat_DD.dat'
out_catalogue='outcat_napa.txt'
Vp=6000 #velocities in m/s !!!!!
Vps=1.73
Vs=Vp/Vps
nref=4
datafile='napa_Vp_'+str(Vp)

dataobj=molquake.molquake_dat(data_path,datafiles,one_sta=True)
locobj=molquake.molquake_loc(Vp,Vs)
references=dataobj.references
nref=num.size(references[:,0])
locobj.location(dataobj.data, references)
catevs=dataobj.read_catalogue(data_path,catalogue_file)
#locobj.catalogue_creation(dataobj.evids, dataobj.origin, nref, out_catalogue, output_path)
if uncertainty:
    multiloc=locobj.uncertainty_estimation(dataobj.data, references, 25, [5500,6500], 1.73)
    locboot=num.array([0.,0.,0])
    for i in range(25):
        locboot=num.vstack((locboot,multiloc[i,:,:]))
    num.save('locboot_real.npy',locboot)
    num.save('locreal.npy',dataobj.references)
    num.save('catevs.npy',catevs)

CVS_E=547327-dataobj.origin[0]
CVS_N=4244262-dataobj.origin[1]

print(num.shape(references),'shape ref')
print(num.shape(locobj.locations),'shape loc')
print(num.shape(dataobj.evids),'shape evids')
nref,mref=num.shape(references)
err_svd=num.sort(num.sqrt((catevs[nref:,0]-locobj.locations[nref:,0])**2+(catevs[nref:,1]-locobj.locations[nref:,1])**2+(catevs[nref:,2]-locobj.locations[nref:,2])**2))
print(err_svd)
n_bins=30
#c1='#FF6347'
#c2='#4682B4'
#c3='#99e699'
msize=30
c1='#4285F4'
c2='#EA4335'
c4='#34A853'
c3='#FBBC05'
nevs=num.size(locobj.locations[:,1])-nref
fig=plt.figure(figsize=(20,20))
ax1=plt.subplot(2,2,1)
#for i in range(nevs):
    #x=[catevs[nref+i,0]+4000,locobj.locations[nref+i,0]+4000]
    #y=[catevs[nref+i,1]+2000,locobj.locations[nref+i,1]+2000]
    #ax1.plot(x,y,color=c2, zorder=2)
ax1.scatter(locobj.locations[nref:,0],locobj.locations[nref:,1], s=50, c=c1, zorder=3)
ax1.scatter(catevs[nref:,0],catevs[nref:,1],s=25, c=c2, zorder=3)
ax1.scatter(dataobj.references[:,0],dataobj.references[:,1], s=50, c=c3, zorder=4)
ax1.set_xlim([-10000,10000])
ax1.set_ylim([-10000,10000])
ax1.set_xticks(range(-10000,12000,2000))
ax1.set_yticks(range(-10000,12000,2000))
ax1.set_xticklabels(range(-10,12,2),fontsize='x-large')
ax1.set_yticklabels(range(-10,12,2),fontsize='x-large')
ax1.set_xlabel('X(km)',fontsize='xx-large')
ax1.set_ylabel('Y(km)',fontsize='xx-large')
ax1.grid('on')
ax1.set_aspect('equal')

ax2=plt.subplot(2,2,2, sharey=ax1)
#for i in range(nevs):
#    x=[catevs[nref+i,2],locobj.locations[nref+i,2]]
#    y=[catevs[nref+i,1],locobj.locations[nref+i,1]]
#    ax2.plot(x,y,color=c2, zorder=2)
ax2.scatter(locobj.locations[nref:,2],locobj.locations[nref:,1],s=50, c=c1, zorder=3)
ax2.scatter(catevs[nref:,2],catevs[nref:,1],s=25, c=c2, zorder=3)
ax2.scatter(dataobj.references[:,2],dataobj.references[:,1],s=50, c=c3, zorder=4)
ax2.set_xlim([0,20000])
ax2.set_ylim([-10000,10000])
ax2.set_xticks(range(0,22000,2000))
ax2.set_yticks(range(-10000,12000,2000))
ax2.set_xticklabels(range(0,22,2),fontsize='x-large')
ax2.set_yticklabels(range(-10,12,2),fontsize='x-large')
ax2.set_xlabel('Z(km)',fontsize='xx-large')
ax2.set_ylabel('Y(km)',fontsize='xx-large')
ax2.set_aspect('equal')
ax2.grid('on')

ax3=plt.subplot(2,2,3, sharex=ax1)
#for i in range(nevs):
#    x=[catevs[nref+i,0],locobj.locations[nref+i,0]]
#    y=[catevs[nref+i,2],locobj.locations[nref+i,2]]
#    ax3.plot(x,y,color=c2, zorder=2)
ax3.scatter(locobj.locations[nref:,0],locobj.locations[nref:,2],s=50, c=c1, zorder=3)
ax3.scatter(catevs[nref:,0],catevs[nref:,2],s=25, c=c2, zorder=3)
ax3.scatter(dataobj.references[:,0],dataobj.references[:,2],s=50, c=c3, zorder=4)
ax3.set_ylim([0,20000])
ax3.set_xlim([-10000,10000])
ax3.set_yticks(range(0,22000,2000))
ax3.set_xticks(range(-10000,12000,2000))
ax3.set_yticklabels(range(0,22,2),fontsize='x-large')
ax3.set_xticklabels(range(-10,12,2),fontsize='x-large')
ax3.set_xlabel('X(km)',fontsize='xx-large')
ax3.set_ylabel('Z(km)',fontsize='xx-large')
ax3.invert_yaxis()
ax3.set_aspect('equal')
ax3.grid('on')

ax4=plt.subplot(2,2,4)
ax4.hist(err_svd, n_bins, density=True, histtype='step', cumulative=True, color=c4)
ax4.set_xlim([0,10000])
ax4.set_ylim([0,1])
ax4.set_xticks(range(0,12000,2000))
ax4.set_yticklabels(['%2.1f'%(i*0.1) for i in range(0,12,2)],fontsize='x-large')
ax4.set_xticklabels(range(0,12,2),fontsize='x-large')
ax4.set_xlabel('Distance from Manual Locations (km)',fontsize='xx-large')
ax4.set_ylabel('Fraction of events',fontsize='xx-large')
ax4.grid('on')

plt.savefig(datafile+'_'+str(nref)+'_small.eps')
plt.show()


# plt.plot(multiloc[:,:,0],multiloc[:,:,1],'ro')
# plt.plot(dataobj.references[:,0],dataobj.references[:,1],'b+')
# plt.xlim([-10000,10000])
# plt.ylim([-15000,15000])
# plt.grid('on')
# plt.show()
#
# plt.plot(multiloc[:,:,0],multiloc[:,:,2],'ro')
# plt.plot(dataobj.references[:,0],dataobj.references[:,2],'b+')
# plt.xlim([-15000,15000])
# plt.ylim([-15000,15000])
# plt.grid('on')
# plt.show()
#
# plt.plot(multiloc[:,:,1],multiloc[:,:,2],'ro')
# plt.plot(dataobj.references[:,1],dataobj.references[:,2],'b+')
# plt.xlim([-15000,15000])
# plt.ylim([-15000,15000])
# plt.grid('on')
# plt.show()

c1='#4285F4'
c2='#EA4335'
c4='#34A853'
c3='#FBBC05'
nevs=num.size(locobj.locations[:,1])-nref
fig=plt.figure(figsize=(20,20))
ax1=plt.subplot(1,1,1)
ax1.scatter(locobj.locations[nref:,0],locobj.locations[nref:,1], s=50, c=c1, zorder=3)
ax1.scatter(catevs[nref:,0],catevs[nref:,1],s=25, c=c2, zorder=3)
ax1.scatter(dataobj.references[:,0],dataobj.references[:,1], s=50, c=c3, zorder=4)
ax1.scatter(CVS_E,CVS_N, marker='v', c=c4,  s=100)
ax1.set_xlim([-15000,15000])
ax1.set_ylim([-15000,15000])
ax1.set_xticks(range(-15000,20000,5000))
ax1.set_yticks(range(-15000,20000,5000))
ax1.set_xticklabels(range(-15,20,5),fontsize='x-large')
ax1.set_yticklabels(range(-15,20,5),fontsize='x-large')
ax1.set_xlabel('X(km)',fontsize='xx-large')
ax1.set_ylabel('Y(km)',fontsize='xx-large')
ax1.grid('on')
ax1.set_aspect('equal')

plt.savefig(datafile+'_'+str(nref)+'.eps')
plt.show()
