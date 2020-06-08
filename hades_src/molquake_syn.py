import molquake
import numpy as num
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plotting(trueloc,molqloc):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(trueloc[:,0],trueloc[:,1],trueloc[:,2],c='b')
    ax.scatter(molqloc[:,0],molqloc[:,1],molqloc[:,2],c='r')
    plt.show()


path='/Users/francesco/Desktop/Single_Station_Location/Synthetic_data'
data_path='/Users/francesco/Desktop/Single_Station_Location/'
filecat='catalogue'
filepicks='picks'
offset=10000
extension=1000
Vp=6500
Vs=Vp/num.sqrt(3)
kv=(Vp*Vs)/(Vp-Vs)
data=num.load(path+'/synth1.data.mol.npy')
references=num.load(path+'/synth1.reference.mol.npy')
references=references[:10,:]
data=num.sqrt(data[:,0]**2+data[:,1]**2+data[:,2]**2)/kv
objloc=molquake.molquake_loc(offset,extension)


Vpmin=6000
Vpmax=7000
locations=objloc.location(data, references, Vp, Vs)
num.save('locations_single.npy',locations)
for i in range(1):
    Vp=Vpmin+(2*num.random.random_sample()-1)*(Vpmax-Vpmin)
    Vs=Vp/num.sqrt(3)
    locs=objloc.location(data, references, Vp, Vs)
    locations=num.vstack((locations,locs))
num.save('locations_bootstrap.npy',locations)
#multiloc=objloc.uncertainty_estimation(data, references, 1, [5.4, 5.6], num.sqrt(3))
#print(num.shape(multiloc))
#data=data_in_sphere(n_points, radius=1000, origin_shift=10000.)
#references=data_in_sphere(n_points, radius=1000, origin_shift=10000.)
data=num.load(path+'/synth1.data.mol.npy')

# for i in range(1):
#     plt.plot(multiloc[i,:,0],multiloc[i,:,1],'ro')
# plt.plot(data[:,0],data[:,1],'b+')
# plt.plot(references[:,0],references[:,1],'g+')
# plt.show()
#
# for i in range(1):
#     plt.plot(multiloc[i,:,0],multiloc[i,:,2],'ro')
# plt.plot(data[:,0],data[:,2],'b+')
# plt.plot(references[:,0],references[:,2],'g+')
# plt.show()
#
# for i in range(1):
#     plt.plot(multiloc[i,:,1],multiloc[i,:,2],'ro')
# plt.plot(data[:,1],data[:,2],'b+')
# plt.plot(references[:,1],references[:,2],'g+')
# plt.show()

plt.plot(locations[:,0],locations[:,1],'ro')
plt.plot(data[:,0],data[:,1],'b+')
plt.plot(references[:,0],references[:,1],'g+')
plt.show()

plt.plot(locations[:,0],locations[:,2],'ro')
plt.plot(data[:,0],data[:,2],'b+')
plt.plot(references[:,0],references[:,2],'g+')
plt.show()

plt.plot(locations[:,1],locations[:,2],'ro')
plt.plot(data[:,1],data[:,2],'b+')
plt.plot(references[:,1],references[:,2],'g+')
plt.show()

nref,mref=num.shape(references)
err_svd=num.sort(num.sqrt((data[:,0]-objloc.locations[nref:,0])**2+(data[:,1]-objloc.locations[nref:,1])**2+(data[:,2]-objloc.locations[nref:,2])**2))

plotting(data,objloc.locations)

n_bins=30
ax=plt.subplot(111)
n1, bins1, patches1 = ax.hist(err_svd, n_bins, density=True, histtype='step', cumulative=True)
plt.savefig('figure_nref_'+str(nref)+'eps')
plt.show()
