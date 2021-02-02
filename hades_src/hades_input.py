import numpy as num
import datetime
import LatLongUTMconversion
import os
import sys

km=1000.

class hades_input:


    def __init__(self,data_path,event_file,station_file):
        events_file=os.path.join(data_path,event_file)
        # event file
        references,refevid,evtsp,reforig,events=hades_input.__read_evfile(event_file)
        self.references=references
        self.refevid=refevid
        self.origin=reforig
        self.data=evtsp
        self.events=events
        # station file
        station_file=os.path.join(data_path,station_file)
        stations=hades_input.__read_stafile(station_file,reforig)
        self.stations=stations


    def __read_evfile(input_file):
        with open(input_file, 'r') as f:
            toks=f.readline().split(';')
            z0,e0,n0=LatLongUTMconversion.LLtoUTM(23, eval(toks[1]), eval(toks[2])) #order lat lon
            refor=(e0,n0,z0)
            refevid=[]
            events=[]
            evtsp={}
            references=[]
            for line in f:
                toks=line.split(';')
                if toks[0][0]=='#':
                    evid=toks[0]
                    evtsp[evid]={}
                    evdate=toks[1][0:10]
                    if toks[0][1]=='R':
                        z,e,n=LatLongUTMconversion.LLtoUTM(23, eval(toks[2]), eval(toks[3])) # order lat lon
                        depth=eval(toks[4])*km
                        refevid.append(evid)
                        references.append([e-e0,n-n0,depth])
                    else:
                        events.append(evid)
                else:
                    sta=toks[0]
                    if toks[1]=='na' or toks[2]=='na':
                        continue
                    tp=datetime.datetime.strptime(evdate+'T'+toks[1], '%Y/%m/%dT%H:%M:%S.%f')
                    ts=datetime.datetime.strptime(evdate+'T'+toks[2], '%Y/%m/%dT%H:%M:%S.%f')
                    tsp=(ts-tp).total_seconds()
                    evtsp[evid][sta]=[tp,ts,tsp]
        references=num.array(references)
        return references,refevid,evtsp,refor,events


    def __read_stafile(input_file,reforig):
        (e0,n0,z0)=reforig
        stations={}
        with open(input_file, 'r') as f:
            for line in f:
                toks=line.split()
                sta=toks[0]
                z,e,n=LatLongUTMconversion.LLtoUTM(23, eval(toks[1]), eval(toks[2])) #order lat lon
                elev=eval(toks[3])*km
                if z==z0:
                    stations[sta]=[e-e0,n-n0,elev]
                else:
                    print('cluster and stations are in different UTM zones')
                    sys.exit()
        return stations


    def distance_calculation(self, Vp, Vs, sta):
        ''' This method calculates the inter-event distance matrix
        for the entire dataset. Input: It requires the P and S wave velocities and
        the station name. For single station mode sta need to be  string with the name of the
        station "STANAME", for multistation sta need to be a list with the stations
        ["STANAME_1", "STANAME_2", ... ,"STANAME_N"].
        '''
        evrefid=self.refevid
        evrefs=self.references
        self.sel_sta=sta
        evids=list((self.data).keys())
        events=self.refevid+self.events
        nevs=len(events)
        nref=len(evrefid)
        distances=num.zeros([nevs,nevs])
        kv=(Vp*Vs)/(Vp-Vs)
        for i in range(nevs-1):
            for j in range(i+1,nevs):
                if (i<nref) and (j<nref):
                   distances[i,j]=num.sqrt((evrefs[i][0]-evrefs[j][0])**2+(evrefs[i][1]-evrefs[j][1])**2+(evrefs[i][2]-evrefs[j][2])**2)
                else:
                   tsp_ev1=self.data[events[i]]
                   tsp_ev2=self.data[events[j]]
                   distances[i,j]=hades_input.__interev_distance(tsp_ev1,tsp_ev2,kv,sta)
                distances[j,i]=distances[i,j]
        self.distances=distances
        self.events=events


    def __interev_distance(tsp_ev1,tsp_ev2,kv,sta):
        if type(sta)==str:
            ie_dist=hades_input.__onesta_interev_distance(tsp_ev1,tsp_ev2,kv,sta)
        elif type(sta)==list and len(sta)==2:
            ie_dist=hades_input.__twosta_interev_distance(tsp_ev1,tsp_ev2,kv,sta)
        elif type(sta)==list and len(sta)>2:
            ie_dist=hades_input.__multi_interev_distance(tsp_ev1,tsp_ev2,kv,sta)
        else:
            print('Error in reading the station list for interevent distance')
            sys.exit()
        return ie_dist


    def __onesta_interev_distance(tsp_ev1,tsp_ev2,kv,sta):
        stalist=list(set(tsp_ev1.keys()) & set(tsp_ev2.keys()))
        if sta in stalist:
            ie_dist=num.abs(tsp_ev1[sta][-1]-tsp_ev2[sta][-1])*kv
            return ie_dist
        else:
            return num.NaN


    def __twosta_interev_distance(tsp_ev1,tsp_ev2,kv,sta):
        sta1=sta[0]
        sta2=sta[1]
        stalist=list(set(tsp_ev1.keys()) & set(tsp_ev2.keys()))
        if sta1 in stalist and sta2 in stalist:
            iedist_sta1=num.abs(tsp_ev1[sta1][-1]-tsp_ev2[sta1][-1])*kv
            iedist_sta2=num.abs(tsp_ev1[sta2][-1]-tsp_ev2[sta2][-1])*kv
            ie_dist=num.sqrt(iedist_sta1**2+iedist_sta2**2)
            return ie_dist
        else:
            return num.NaN


    def __multi_interev_distance(tsp_ev1,tsp_ev2,kv,sta):

        stalist=list((set(tsp_ev1.keys()) & set(tsp_ev2.keys())) & set(sta))

        nsta=len(stalist)

        if nsta>=2:
            ie_dist=[];
            for i in range(nsta-1):
                sta1=stalist[i]
                iedist_sta1=num.abs(tsp_ev1[sta1][-1]-tsp_ev2[sta1][-1])*kv
                for j in range(i+1,nsta):
                    sta2=stalist[j]
                    iedist_sta2=num.abs(tsp_ev1[sta2][-1]-tsp_ev2[sta2][-1])*kv
                    ie_dist.append(num.sqrt(iedist_sta1**2+iedist_sta2**2))
            ie_dist=num.mean(num.array(ie_dist))

            return ie_dist

        else:

            return num.NaN


# if __name__ == "__main__":
#     hobj1=hades_input('./','Sequenza_KT_2020.csv','stations.txt')
#     hobj2=hades_input('./','Sequenza_KT_2020.csv','stations.txt')
#     hobj1.distance_calculation(5000,3000,['FRAE','SEME'],)
#     hobj2.distance_calculation(5000,3000,['SEME','FRAE','GIZE','MUSE'])
#     for i in range(len(hobj1.events)):
#         for j in range(i+1,len(hobj1.events)):
#             print(hobj1.events[i],hobj1.events[j],'2 sta',hobj1.distances[i,j],'4 sta',hobj2.distances[i,j])
