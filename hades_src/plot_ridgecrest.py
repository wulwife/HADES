import numpy as num
import matplotlib.pyplot as plt
import scipy.stats as st

xmin=-25000; xmax=25000
ymin=-25000; ymax=25000
#edata=num.load('locations_bootstrap.npy')
#tdata=num.load('synth1.data.mol.npy')
edata=num.load('locboot_real_ridgecrest.npy')
tdata=num.load('catevs_ridgecrest.npy')
ksmoothing=0.5
xx, yy = num.mgrid[xmin:xmax:100j, ymin:ymax:100j]
xe=edata[:,0]+4000; ye=edata[:,1]+2000
xt=tdata[:,0]+4000; yt=tdata[:,1]+2000
positions = num.vstack([xx.ravel(), yy.ravel()])
values = num.vstack([xe, ye])
kernel = st.gaussian_kde(values, ksmoothing)
f = num.reshape(kernel(positions).T, xx.shape)


fig = plt.figure(figsize=(10,10))
ax = fig.gca()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
ax.imshow(num.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
cset = ax.contour(xx, yy, f, colors='k')
#ax.scatter(xe,ye,c='r',s=10)
ax.scatter(xt,yt,c='w',s=10)
#ax.clabel(cset, inline=1, fontsize=10)
plt.title('2D Gaussian Kernel density estimation')
plt.axis('equal')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks(range(xmin,xmax+5000,5000))
ax.set_yticks(range(ymin,ymax+5000,5000))
for i in range(xmin//1000,(xmax+5000)//1000):
    print(i)
ax.set_xticklabels(range(xmin//1000,(xmax+5000)//1000,5),fontsize='x-large')
ax.set_yticklabels(range(ymin//1000,(ymax+5000)//1000,5),fontsize='x-large')
ax.set_xlabel('X(km)',fontsize='xx-large')
ax.set_ylabel('Y(km)',fontsize='xx-large')
plt.grid('on')
plt.savefig('KDF_XY_ridge.eps')
plt.show()


xmin=-25000; xmax=25000
ymin=0; ymax=50000
xx, yy = num.mgrid[xmin:xmax:100j, ymin:ymax:100j]
xe=edata[:,0]+4000; ye=edata[:,2]
xt=tdata[:,0]+4000; yt=tdata[:,2]
positions = num.vstack([xx.ravel(), yy.ravel()])
values = num.vstack([xe, ye])
kernel = st.gaussian_kde(values, ksmoothing)
f = num.reshape(kernel(positions).T, xx.shape)

fig = plt.figure(figsize=(10,10))
ax = fig.gca()
cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
ax.imshow(num.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
cset = ax.contour(xx, yy, f, colors='k')
#ax.scatter(xe,ye,c='r',s=10)
ax.scatter(xt,yt,c='w',s=10)
#ax.clabel(cset, inline=1, fontsize=10)
ax.set_xlabel('X(km)',fontsize='xx-large')
ax.set_ylabel('Z(km)',fontsize='xx-large')
plt.title('2D Gaussian Kernel density estimation')
plt.axis('equal')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks(range(xmin,xmax+5000,5000))
ax.set_yticks(range(ymin,ymax+5000,5000))
ax.set_xticklabels(range(xmin//1000,(xmax+5000)//1000,5),fontsize='x-large')
ax.set_yticklabels(range(ymin//1000,(ymax+5000)//1000,5),fontsize='x-large')
ax.invert_yaxis()
plt.grid('on')
plt.savefig('KDF_XZ_ridge.eps')
plt.show()

xmin=-25000; xmax=25000
ymin=0; ymax=50000
xx, yy = num.mgrid[xmin:xmax:100j, ymin:ymax:100j]
xe=edata[:,1]+2000; ye=edata[:,2]
xt=tdata[:,1]+2000; yt=tdata[:,2]
positions = num.vstack([xx.ravel(), yy.ravel()])
values = num.vstack([xe, ye])
kernel = st.gaussian_kde(values, ksmoothing)
f = num.reshape(kernel(positions).T, xx.shape)

fig = plt.figure(figsize=(10,10))
ax = fig.gca()
cfset = ax.contourf(xx, yy, f, cmap='coolwarm')
ax.imshow(num.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
cset = ax.contour(xx, yy, f, colors='k')
#ax.scatter(xe,ye,c='r',s=10)
ax.scatter(xt,yt,c='w',s=10)
#ax.clabel(cset, inline=1, fontsize=10)
ax.set_xlabel('Y(km)',fontsize='xx-large')
ax.set_ylabel('Z(km)',fontsize='xx-large')
plt.title('2D Gaussian Kernel density estimation')
plt.axis('equal')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks(range(xmin,xmax+5000,5000))
ax.set_yticks(range(ymin,ymax+5000,5000))
ax.set_xticklabels(range(xmin//1000,(xmax+5000)//1000,5),fontsize='x-large')
ax.set_yticklabels(range(ymin//1000,(ymax+5000)//1000,5),fontsize='x-large')
ax.invert_yaxis()
plt.grid('on')
plt.savefig('KDF_YZ_ridge.eps')
plt.show()
