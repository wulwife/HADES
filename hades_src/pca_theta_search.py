import numpy as num
import LatLongUTMconversion
import matplotlib.pyplot as plt
from hades_input import hades_input


def read_catalogue(filename,hobj):
    e0,n0,z0=hobj.origin[0],hobj.origin[1],hobj.origin[2]
    x_evs=[]; y_evs=[]; z_evs=[]; evtsps={}
    #Id Id Lat Lon Depth Error Station(s) Date Tp Ts-Tp Station(s) Date Tp Ts-Tp
    f=open(filename, 'r')
    f.readline()

    for line in f:
        toks=line.split()
        lat_ev, lon_ev, dep_ev = eval(toks[1]), eval(toks[2]), eval(toks[3])
        sta_name1, tsp_sta1 = toks[4], eval(toks[5])
        sta_name2, tsp_sta2 = toks[8], eval(toks[9])
        z,e,n=LatLongUTMconversion.LLtoUTM(23, lat_ev, lon_ev)
        x_evs.append(e-e0)
        y_evs.append(n-n0)
        z_evs.append(dep_ev*1000)
        if (sta_name1 not in evtsps.keys()) and (sta_name2 not in evtsps.keys()):
            evtsps[sta_name1]=[]
            evtsps[sta_name2]=[]
        evtsps[sta_name1].append(tsp_sta1)
        evtsps[sta_name2].append(tsp_sta2)
    f.close()

    x_evs=num.array(x_evs)
    y_evs=num.array(y_evs)
    z_evs=num.array(z_evs)

    return x_evs, y_evs, z_evs, evtsps



def pca_theta_search(filename,hobj,angle_increment):
    x_evs, y_evs, z_evs, evtsps = read_catalogue(filename,hobj)
    c_evs=x_evs+1j*y_evs
    rect_max=0
    rects=[]
    angle_increment=int(angle_increment)
    thetas=num.arange(0,360+angle_increment,angle_increment)
    for theta in thetas:
        theta_r=theta*(num.pi/180)
        c_rot=c_evs*num.exp(-1j*theta_r)
        rect=1
        data={}
        stas=list((hobj.stations).keys())
        sta1=stas[0]
        sta2=stas[1]
        for sta in (hobj.stations).keys():
            X=num.zeros([num.size(x_evs),2])
            data[sta]=num.zeros([num.size(x_evs),2])
            dx=(c_rot.real-hobj.stations[sta][0])
            dy=(c_rot.imag-hobj.stations[sta][1])
            dz=(z_evs-hobj.stations[sta][2])
            tsp=num.array(evtsps[sta])
            #tsp=tsp_obs-num.mean(tsp_obs)
            dist=num.sqrt(dx**2+dy**2+dz**2)
            ir_dist=num.argsort(dist)
            X[:,0]=dist[ir_dist]
            X[:,1]=tsp[ir_dist]
            data[sta][:,0]=dist[ir_dist]
            data[sta][:,1]=tsp[ir_dist]
            M=num.mean(X.T, axis=1)
            C=X-M
            V=num.cov(C.T)
            values, vectors = num.linalg.eig(V)
            values=values/num.max(values)
            rect=rect*(num.max(values)/num.min(values))
        rects.append(rect)
        plt.title(str(theta)+' rect value : '+str(rect))
        plt.plot(data[sta1][:,0],data[sta1][:,1],'or')
        plt.plot(data[sta2][:,0],data[sta2][:,1],'ob')
        plt.show()
        if rect>rect_max:
            rect_max=rect
            theta_best=theta_r
    rects=num.array(rects)
    crot=c_evs*num.exp(-1j*theta_best)
    x_evs=crot.real
    y_evs=crot.imag
    return x_evs, y_evs, rects




#lat,lon=LatLongUTMconversion.UTMtoLL(23, self.locations[i,1]+hobj.origin[1], self.locations[i,0]+hobj.origin[0],hobj.origin[2])



cat_filename='hades_castor_abs.txt'
data_path='./'
input_file='castor_refstich.dat'#_4_sta.csv'#'ridgecrest_gji.dat'
sta_file='stations.txt'#'stations_ridge.txt'
out_file='hades_castor'#'hades_ridgecrest'
stations=['ALCN','ALCX']
hobj=hades_input(data_path,input_file,sta_file)

angle_increment=10
x_evs, y_evs, rects=pca_theta_search(cat_filename,hobj,angle_increment)
thetas=num.arange(0,360+angle_increment,angle_increment)
rects=rects/num.max(rects)
print(num.size(rects),num.size(thetas))
plt.plot(thetas,rects)
plt.show()
