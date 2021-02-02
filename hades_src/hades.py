import numpy as num
from hades_input import hades_input
from hades_location import hades_location
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import LatLongUTMconversion

def tplot(ifig):
    fig=plt.figure(figsize=(10.0,15.0))
    ax=plt.subplot(111)
    lat_evs=[]; lon_evs=[]; frae_evs=[]; seme_evs=[]
    nrm=Normalize(vmin=3., vmax=7.0)
    cmap = plt.cm.get_cmap('RdBu', 20)
    z0,e0,n0=LatLongUTMconversion.LLtoUTM(23, 39.090667, 17.100000)
    f=open('test'+str(ifig)+'.txt', 'r')
    f.readline()
    for line in f:
        toks=line.split()
        lat_ev, lon_ev, seme, frae=eval(toks[1]), eval(toks[2]), eval(toks[-3]), eval(toks[-1])
        z,e,n=LatLongUTMconversion.LLtoUTM(23, lat_ev, lon_ev)
        lat_evs.append(n-n0)
        lon_evs.append(e-e0)
        seme_evs.append(seme)
        frae_evs.append(frae)
    f.close()

    xst=num.zeros(4)
    yst=num.zeros(4)
    yst[0]=39.168111; xst[0]= 17.097250
    yst[1]=39.057294; xst[1]= 17.075408
    yst[2]=39.024306; xst[2]= 17.199000
    yst[3]=38.994306; xst[3]= 17.150306

    for i in range(4):
        z,es,ns=LatLongUTMconversion.LLtoUTM(23,yst[i],xst[i])
        ax.scatter(es-e0, ns-n0, c='#eed54f', marker='v', s=200, zorder=3, linewidth=0.5)

    ax.scatter(lon_evs, lat_evs, c=seme_evs, cmap=cmap)
    plt.grid('on')
    plt.savefig('figure'+str(ifig)+'.eps')
    return None


Vp=5300
Vs=Vp/1.73
hobj=hades_input('./','Sequenza_KT_2020_mod.csv','stations.txt')
hobj.distance_calculation(Vp,Vs,['SEME','FRAE'])
hloc=hades_location(hobj,'./')
hloc.location('test'+str(100)+'.txt')
tplot(100)
