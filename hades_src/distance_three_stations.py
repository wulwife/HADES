import molquake
import numpy as num
import matplotlib.pyplot as plt

def interev_distance(stations,tsp_ev1,tsP_ev2):

    stalist=stalist(stations.keys())
    R={}
    for sta in stations.keys():
        r1=distance_from_tsp(tsp_ev1[sta])**2
        r2=distance_from_tsp(tsp_ev2[sta])**2
        R[sta]=[r1,r2]

    H=[]; S=[]
    for sta1 in stations.keys():
        for sta2 in stations.keys():
            if sta1!=sta2:
                H.append((R[sta1][0]-R[sta2][0])-(R[sta1][1]-R[sta2][1])/2.)
                S.append([(stations[sta2][0]-stations[sta1][0]),
                          (stations[sta2][1]-stations[sta1][1]),
                          (stations[sta2][2]-stations[sta1][2])])

    H=num.array(H)
    S=num.array(S)

    GT=num.dot(S.T,S)
    Ginv=num.linalg.inv(GT)
    GG=num.dot(Ginv,S.T)
    m=num.dot(GG,H)
    ie_dist=num.sqrt(num.sum(m**2))

    return ie_dist


def distance_from_tsp(tsp,Vp=4700,Vps=1.73):
    Vs=Vp/Vps
    k=((Vp*Vs)/(Vp-Vs))
    distance=k*tsp
    return distance
