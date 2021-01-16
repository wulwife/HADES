import molquake
import numpy as num
import matplotlib.pyplot as plt

def interev_distance(tsp_ev1,tsp_ev2):

    stalist=list(set(tsp_ev1.keys()) & set(tsp_ev2.keys()))
    nsta=len(stalist)

    R={}
    for sta in stalist:
        r1=distance_from_tsp(tsp_ev1[sta])**2
        r2=distance_from_tsp(tsp_ev2[sta])**2
        print(sta,r1,r2)
        R[sta]=[r1,r2]

    H=[]; S=[]
    for i in range(nsta-1):
        for j in range(i+1,nsta):
            sta1=stalist[i]
            sta2=stalist[j]
            H.append(((R[sta1][0]-R[sta1][1])-(R[sta2][0]-R[sta2][1]))/2.)
            S.append([(stations[sta2][0]-stations[sta1][0]),
                      (stations[sta2][1]-stations[sta1][1]),
                      (stations[sta2][2]-stations[sta1][2])])

    H=num.array(H)
    S=num.array(S)

    GT=num.dot(S.T,S)
    Ginv=num.linalg.inv(GT)
    GG=num.dot(Ginv,S.T)
    #GG=num.linalg.pinv(S)
    m=num.dot(GG,H)
    ie_dist=num.sqrt(num.sum(m**2))

    return ie_dist


def distance_from_tsp(tsp,Vp=4700,Vps=1.73):
    Vs=Vp/Vps
    k=((Vp*Vs)/(Vp-Vs))
    distance=k*tsp
    return distance

if __name__ == "__main__":
    stations={'sta1':num.array([0,0,0]),
    'sta2':num.array([0,10,0]),
    'sta3':num.array([10,0,0]),
    'sta4':num.array([0,0,10])}

    Vp=4700
    Vps=1.73
    Vs=Vp/Vps
    k=((Vp*Vs)/(Vp-Vs))

    ev1=num.array([50,50,1050])
    ev2=num.array([25,25,2050])

    tsp_ev1={}
    tsp_ev2={}

    for sta in stations.keys():
        tsp_ev1[sta]=num.sqrt(num.sum((stations[sta]-ev1)**2))/k
        tsp_ev2[sta]=num.sqrt(num.sum((stations[sta]-ev2)**2))/k

    obs_dist=interev_distance(tsp_ev1,tsp_ev2)
    cal_dist=num.sqrt(num.sum((ev2-ev1)**2))
    print('results',cal_dist,obs_dist)
