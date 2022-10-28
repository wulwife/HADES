import numpy as num
import datetime
import latlon2cart
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
            latref,lonref=eval(toks[1]),eval(toks[2])
            orig=latlon2cart.Coordinates(latref,lonref,0)
            try:
                depthref=eval(toks[3])*km
            except:
                depthref=0.
            refor=(latref,lonref,depthref)
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
                        e,n,z = orig.geo2cart(eval(toks[2]), eval(toks[3]),0)
                        depth=eval(toks[4])*km
                        refevid.append(evid)
                        references.append([e,n,depth])
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


    def __read_stafile(input_file,refor):
        (latref,lonref,depthref)=refor
        orig=latlon2cart.Coordinates(latref,lonref,0)
        stations={}
        with open(input_file, 'r') as f:
            for line in f:
                toks=line.split()
                sta=toks[0]
                e,n,z = orig.geo2cart(eval(toks[1]), eval(toks[2]),0)
                elev=eval(toks[3])*km
                stations[sta]=[e,n,elev]
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
                   distances[i,j]=hades_input.__interev_distance(tsp_ev1,tsp_ev2,kv,sta,self.stations)
                distances[j,i]=distances[i,j]
        self.distances=distances
        self.events=events
        self.vp=Vp
        self.vs=Vs


    def __interev_distance(tsp_ev1,tsp_ev2,kv,sta,stations):
        if (type(sta)==str) and sta!='ALL':
            ie_dist=hades_input.__onesta_interev_distance(tsp_ev1,tsp_ev2,kv,sta)
        elif type(sta)==list and len(sta)==2:
            ie_dist=hades_input.__twosta_interev_distance(tsp_ev1,tsp_ev2,kv,sta)
        elif (type(sta)==list and len(sta)>2) or sta=='ALL':
            ie_dist=hades_input.__multi_interev_distance(tsp_ev1,tsp_ev2,kv,sta,stations)
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


    def __multi_interev_distance(tsp_ev1,tsp_ev2,kv,sta,stations):

        if sta=='ALL':
            stalist=list(set(tsp_ev1.keys()) & set(tsp_ev2.keys()))
        else:
            stalist=list(set(tsp_ev1.keys()) & set(tsp_ev2.keys()) & set(sta))

        nsta=len(stalist)

        # R={}
        # for ista in stalist:
        #      r1=(num.abs(tsp_ev1[ista][-1])*kv)**2
        #      r2=(num.abs(tsp_ev2[ista][-1])*kv)**2
        #      R[ista]=[r1,r2]


        H=[]; S=[]
        for i in range(nsta-1):
            sta1=stalist[i]
            for j in range(i+1,nsta):
                sta2=stalist[j]
                dsta1=num.abs((tsp_ev1[sta1][-1])**2-num.abs(tsp_ev2[sta1][-1])**2)
                dsta2=num.abs((tsp_ev1[sta2][-1])**2-num.abs(tsp_ev2[sta2][-1])**2)
                H.append(((dsta1-dsta2)*kv**2)/2.)
                #H.append(((R[sta1][0]-R[sta1][1])-(R[sta2][0]-R[sta2][1]))/2.)
                S.append([(stations[sta2][0]-stations[sta1][0]),
                        (stations[sta2][1]-stations[sta1][1]),
                        (stations[sta2][2]-stations[sta1][2])])

        H=num.array(H)
        S=num.array(S)

        GT=num.dot(S.T,S)
        Ginv=num.linalg.inv(GT)
        GG=num.dot(Ginv,S.T)
        GG=num.linalg.pinv(S) #pasudo inverse
        m=num.dot(GG,H)
        ie_dist=num.sqrt(num.sum(m**2))

        return ie_dist


    def relative_frame(self,Vp,Vs,sta,y_ref=-1,z_ref=-1,fixed_depth=0):

        kv=(Vp*Vs)/(Vp-Vs)
        if len(self.refevid)>4:
            events=self.refevid[0:4]
        else:
            events=self.refevid

        d=num.zeros([4,4])
        for i in range(3):
            tsp_ev1=self.data[events[i]]
            for j in range(i+1,4):
                tsp_ev2=self.data[events[j]]
                d[i,j]=hades_input.__interev_distance(tsp_ev1,tsp_ev2,kv,sta,self.stations)
                d[j,i]=d[i,j]

        references=num.zeros([4,3])

        references[0,0]=0.
        references[0,1]=0.
        references[0,2]=0.

        references[1,0]=d[0,1]
        references[1,1]=0.
        references[1,2]=0.

        references[2,0]=(d[0,2]**2-d[1,2]**2)/(2*references[1,0])+(references[1,0]/2)
        references[2,1]=y_ref*num.sqrt(d[0,2]**2-references[2,0]**2)
        references[2,2]=0.

        if fixed_depth:
            references[3,0]=(d[0,3]**2-d[1,3]**2)/(2*references[1,0])+(references[1,0]/2)
            references[3,1]=(d[1,3]**2-d[2,3]**2-(references[3,0]-references[1,0])**2+(references[3,0]-references[2,0])**2)/(2*references[2,1])+(references[2,1]/2)
            references[3,2]=fixed_depth*1000.
        else:
            dmax=num.max(d[:,3])
            xmax=num.max(num.abs(references[:,0]))+dmax
            ymax=num.max(num.abs(references[:,1]))+dmax
            xax=-xmax+((num.arange(100)*0.01)*2*xmax)
            yax=-ymax+((num.arange(100)*0.01)*2*ymax)
            zax=((num.arange(50)*0.02)*dmax)
            errmin=1E10
            for x in xax:
                for y in yax:
                    for z in zax:
                        error=num.sum(num.sqrt((references[0:3,0]-x)**2+(references[0:3,1]-y)**2+(references[0:3,2]-z)**2))
                        if error<errmin:
                            errmin=error
                            references[3,0]=x
                            references[3,1]=y
                            references[3,2]=z

        self.rel_references=references
        self.refevid=events

    def relative_frame_from_file(self):

        if len(self.refevid)>4:
            events=self.refevid[0:4]
            references=self.references[0:4,:]
        else:
            events=self.refevid
            references=self.references

        references[:,2]=references[:,2]-self.origin[-1]
        self.rel_references=references
        self.refevid=events
