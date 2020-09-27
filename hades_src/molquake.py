import numpy as num
import datetime
import LatLongUTMconversion
import os
import sys

km=1000.

class molquake_dat:

    def __init__(self,data_path,datafiles,one_sta=False):
        if type(datafiles)==list and len(datafiles)==2:
          ref_file=os.path.join(data_path,datafiles[0])
          pick_file=os.path.join(data_path,datafiles[1])
          #try:
          if one_sta:
               evtp,evtsp=molquake_dat.read_picksdata(pick_file)
               refevs, refloc, reftsp, reforig=molquake_dat.read_references(ref_file)
          else:
               evtp,evtsp=molquake_dat.read_picksdata_2stations(pick_file)
               refevs, refloc, reftsp, reforig=molquake_dat.read_references_2stations(ref_file)
          #except:
            #print('Catalogue and/or Pick file file not extisitng or not correctly read!')
            #sys.exit()
          nrefshp=num.shape(reftsp)
          if num.size(nrefshp)>1:
                nreftsp,mreftsp=nrefshp
          refs=['REFERENCE' for i in range(nreftsp)]
          self.refcat=refevs
          self.references=refloc
          self.origin=reforig
          self.evids=refs+evtp
          if one_sta:
             self.data=num.hstack((reftsp,evtsp))
          else:
             self.data=num.vstack((reftsp,evtsp))
             print(num.shape(self.data))
        else:
          self.load_data(data_path,datafiles)

    def read_references(ref_file):
        '''read the file containing the reference events,
           it is a text file with five columns as follow:
           origin_time lat lon depth(km) tP tS
           the format for origin_time, tP and tS is the following: yyyymmddThhmmss.ms
           the first row of the file contains information of the origin of reference system
           and it is denoted with the word ORIGIN, followed by latitute and longitude values
           EXAMPLE:
           ORIGIN 38.226782 -122.331413
           20140824T102044.07 38.21517 -122.31233 11.12 20140824T102048.270000 20140824T102051.240000
           20140824T102444.24 38.25983 -122.33733 10.34 20140824T102447.510000 20140824T102449.840000
           20140824T102725.97 38.23917 -122.33783  6.77 20140824T102729.190000 20140824T102731.570000
        '''
        with open(ref_file, 'r') as f:
            refevs=[]
            reftsp=[]
            reflocs=[]
            z0,e0,n0=[0,0,0]
            reforig=[0,0,0]
            for line in f:
                toks=line.split()
                if toks[0]=='ORIGIN' or toks[0]=='origin':
                   z0,e0,n0=LatLongUTMconversion.LLtoUTM(23, eval(toks[1]), eval(toks[2]))
                   reforig=(e0,n0,z0)
                else:
                   otime, lat, lon, depth=toks[0], eval(toks[1]), eval(toks[2]), eval(toks[3])*km,
                   tp=datetime.datetime.strptime(toks[4], '%Y%m%dT%H%M%S.%f')
                   ts=datetime.datetime.strptime(toks[5], '%Y%m%dT%H%M%S.%f')
                   tsp=(ts-tp).total_seconds()
                   z,e,n=LatLongUTMconversion.LLtoUTM(23, lat, lon)
                   reflocs.append([e-reforig[0],n-reforig[1],depth])
                   reftsp.append(tsp)
                   refevs.append([otime,lat,lon,e-reforig[0],n-reforig[1],depth])
        refloc=num.array(reflocs)
        reftsp=num.array(reftsp)
        return refevs, refloc, reftsp, reforig

    def read_references_2stations(ref_file):
        '''read the file containing the reference events,
           it is a text file with five columns as follow:
           origin_time lat lon depth(km) tP tS
           the format for origin_time, tP and tS is the following: yyyymmddThhmmss.ms
           the first row of the file contains information of the origin of reference system
           and it is denoted with the word ORIGIN, followed by latitute and longitude values
           EXAMPLE:
           ORIGIN 38.226782 -122.331413
           20140824T102044.07 38.21517 -122.31233 11.12 20140824T102048.270000 20140824T102051.240000
           20140824T102444.24 38.25983 -122.33733 10.34 20140824T102447.510000 20140824T102449.840000
           20140824T102725.97 38.23917 -122.33783  6.77 20140824T102729.190000 20140824T102731.570000
        '''
        with open(ref_file, 'r') as f:
            refevs=[]
            reftsp=[]
            reflocs=[]
            z0,e0,n0=[0,0,0]
            reforig=[0,0,0]
            for line in f:
                toks=line.split()
                if toks[0]=='ORIGIN' or toks[0]=='origin':
                   z0,e0,n0=LatLongUTMconversion.LLtoUTM(23, eval(toks[1]), eval(toks[2]))
                   reforig=(e0,n0,z0)
                else:
                   otime, lat, lon, depth=toks[0], eval(toks[1]), eval(toks[2]), eval(toks[3])*km,
                   tp1=datetime.datetime.strptime(toks[4], '%Y%m%dT%H%M%S.%f')
                   ts1=datetime.datetime.strptime(toks[5], '%Y%m%dT%H%M%S.%f')
                   tsp1=(ts1-tp1).total_seconds()
                   tp2=datetime.datetime.strptime(toks[6], '%Y%m%dT%H%M%S.%f')
                   ts2=datetime.datetime.strptime(toks[7], '%Y%m%dT%H%M%S.%f')
                   tsp2=(ts2-tp2).total_seconds()
                   z,e,n=LatLongUTMconversion.LLtoUTM(23, lat, lon)
                   reflocs.append([e-reforig[0],n-reforig[1],depth])
                   reftsp.append([tsp1,tsp2])
                   refevs.append([otime,lat,lon,e-reforig[0],n-reforig[1],depth])
        refloc=num.array(reflocs)
        reftsp=num.array(reftsp)
        return refevs, refloc, reftsp, reforig

    def read_picksdata(pick_file):
        '''read the file containing the picks of the events to locate,
           it is a text file with two column,
           the first one for tP and the second one for tS
           the format for tP andtS is the following: yyyymmddThhmmss.ms
           EXAMPLE:
           20140824T105750.650000 20140824T105753.670000
           20140824T110441.320000 20140824T110444.130000
           20140824T110733.830000 20140824T110736.880000
           20140824T111552.660000 20140824T111555.430000
           ...
        '''
        with open(pick_file, 'r') as f:
            evtp=[]
            evtsp=[]
            for line in f:
                toks=line.split()
                if toks[0]=='#':
                    continue
                tp=datetime.datetime.strptime(toks[0], '%Y%m%dT%H%M%S.%f')
                ts=datetime.datetime.strptime(toks[1], '%Y%m%dT%H%M%S.%f')
                tsp=(ts-tp).total_seconds()
                evtp.append(tp)
                evtsp.append(tsp)
        evtsp=num.array(evtsp)
        return evtp, evtsp

    def read_picksdata_2stations(pick_file):
        '''read the file containing the picks of the events to locate,
           it is a text file with two column,
           the first one for tP and the second one for tS
           the format for tP andtS is the following: yyyymmddThhmmss.ms
           EXAMPLE:
           20140824T105750.650000 20140824T105753.670000 20140824T105750.650000 20140824T105753.670000
           20140824T110441.320000 20140824T110444.130000 20140824T110441.320000 20140824T110444.130000
           20140824T110733.830000 20140824T110736.880000 20140824T110733.830000 20140824T110736.880000
           20140824T111552.660000 20140824T111555.430000 20140824T111552.660000 20140824T111555.430000
           ...
        '''
        with open(pick_file, 'r') as f:
            evtp=[]
            evtsp=[]
            for line in f:
                toks=line.split()
                if toks[0]=='#':
                    continue
                tp1=datetime.datetime.strptime(toks[0], '%Y%m%dT%H%M%S.%f')
                ts1=datetime.datetime.strptime(toks[1], '%Y%m%dT%H%M%S.%f')
                tsp1=(ts1-tp1).total_seconds()
                tp2=datetime.datetime.strptime(toks[2], '%Y%m%dT%H%M%S.%f')
                ts2=datetime.datetime.strptime(toks[3], '%Y%m%dT%H%M%S.%f')
                tsp2=(ts2-tp2).total_seconds()
                evtp.append(num.min([tp1, tp2]))
                evtsp.append([tsp1, tsp2])
        evtsp=num.array(evtsp)
        return evtp, evtsp

    def read_catalogue(self,data_path,catalogue_file):
        cat_file=os.path.join(data_path,catalogue_file)
        with open(cat_file, 'r') as f:
            catlocs=[]
            for line in f:
                toks=line.split()
                otime, lat, lon, depth=toks[0], eval(toks[1]), eval(toks[2]), eval(toks[5])*km
                z,e,n=LatLongUTMconversion.LLtoUTM(23, lat, lon)
                catlocs.append([e-self.origin[0],n-self.origin[1],depth])
        catloc=num.array(catlocs)
        return catloc

    def store_data(self,out_path,filename):
        f_out=os.path.join(out_path,filename)
        try:
          num.savez(f_out+'.molquake_data', data=self.data, evids=self.evids, references=self.references, refcat=self.refcat, orig=self.origin)
        except:
          print('Wrong path!!! I will use the current directory to store the files!!')
          num.savez(filename+'.molquake_data', data=self.data, evids=self.evids, references=self.references, refcat=self.refcat, orig=self.origin)

    def load_data(self,data_path,filename):
        f_in=os.path.join(data_path,filename)
        try:
          molquake_data=num.load(f_in)
          self.references=molquake_data['references']
          self.data=molquake_data['data']
          self.refcat=molquake_data['refcat']
          self.evids=molquake_data['refids']
          self.origin=molquake_data['orig']
        except:
          print('file: '+f_in+' not found!')
          sys.exit()



class molquake_loc:

    def __init__(self, Vp, Vs):
        if Vp > Vs:
           self.Vp=Vp
           self.Vs=Vs
        else:
           print('Vp must be larger than Vs')
           sys.exit()

    def dgs_locator(event, references, distances):
        '''event is the id of the event you want locate
        references is an object array of the form ['eventid',x,y,z]
        this method returns the event location and the updated the reference locations
        that include the new event located'''

        n_ref,m_ref=num.shape(references)

        X=num.array([references[:,0],references[:,1],references[:,2]]).T
        XC=num.mean(X,axis=0)
        D=num.zeros([n_ref,n_ref])
        for i in range(n_ref):
            for j in range(n_ref):
                Xi=references[i,0]; Xj=references[j,0]; #need to be optimized
                Yi=references[i,1]; Yj=references[j,1];
                Zi=references[i,2]; Zj=references[j,2];
                dio=distances[i,event]; djo=distances[j,event];
                dij=num.sqrt((Xi-Xj)**2+(Yi-Yj)**2+(Zi-Zj)**2)
                D[i,j]=(dio**2+djo**2-dij**2)/2.
                D[j,i]=D[i,j]

        U,S,V=num.linalg.svd(D, full_matrices=True)
        S=num.diag(S)
        Y=num.dot(U[:,0:3],num.sqrt(S[0:3,0:3]))
        XC=num.mean(X,axis=0)
        YC=num.mean(Y,axis=0)

        XTY=num.dot((X-XC).T,(Y-YC))
        U1,S1,V1=num.linalg.svd(XTY, full_matrices=True)
        Q=num.dot(U1,V1)
        Y=num.dot(Q,(Y-YC).T)
        #X=X-XC
        Xfin=num.dot(-Q,YC)+XC
        evloc=num.array([Xfin[0], Xfin[1], Xfin[2]])
        references=num.vstack((references,evloc))
        return references

    def distance_calculation(data, references, Vp, Vs):
        nref,mref=num.shape(references)
        nevs=num.size(data)
        distances=num.zeros([nevs,nevs])
        kv=(Vp*Vs)/(Vp-Vs)
        for i in range(nevs-1):
            for j in range(i+1,nevs):
                if (i<nref) and (j<nref):
                   distances[i,j]=num.sqrt((references[i,0]-references[j,0])**2+(references[i,1]-references[j,1])**2+(references[i,2]-references[j,2])**2)
                else:
                   distances[i,j]=(num.abs(data[i]-data[j])*kv)
                distances[j,i]=distances[i,j]
        return distances

    def distance_calculation_2stations(data, references, Vp, Vs):
        nref,mref=num.shape(references)
        nevs,mevs=num.shape(data)
        distances=num.zeros([nevs,nevs])
        kv=(Vp*Vs)/(Vp-Vs)
        for i in range(nevs-1):
            for j in range(i+1,nevs):
                if (i<nref) and (j<nref):
                   distances[i,j]=num.sqrt((references[i,0]-references[j,0])**2+(references[i,1]-references[j,1])**2+(references[i,2]-references[j,2])**2)
                else:
                   dist_sta1=num.abs(data[i,0]-data[j,0])*kv
                   dist_sta2=num.abs(data[i,1]-data[j,1])*kv
                   distances[i,j]=num.sqrt(dist_sta1**2+dist_sta2**2)
                distances[j,i]=distances[i,j]
        return distances

    def location_2stations(self, data, references):
        distances=molquake_loc.distance_calculation_2stations(data, references, self.Vp, self.Vs)
        nref,mref=num.shape(references)
        nevs,mevs=num.shape(distances)
        for evid in range(nref,nevs):
            references=molquake_loc.dgs_locator(evid, references, distances)
        self.locations=references
        return references

    def location(self, data, references):
        distances=molquake_loc.distance_calculation(data, references, self.Vp, self.Vs)
        nref,mref=num.shape(references)
        nevs,mevs=num.shape(distances)
        for evid in range(nref,nevs):
            references=molquake_loc.dgs_locator(evid, references, distances)
        self.locations=references
        return references

    def uncertainty_estimation(self, data, references, niters, Vp, Vps=num.sqrt(3)): #not yet complete
        nref,mref=num.shape(references)
        multiloc=[]
        if len(Vp)!=2:
            print('Give a list with only 2 elements for Vp!')
        Vpmin=num.min(Vp); Vpmax=num.max(Vp)
        for niter in range(niters): #this can be parallelized
            Vp=Vpmin+(2*num.random.random_sample()-1)*(Vpmax-Vpmin)
            Vs=Vp/Vps
            distances=molquake_loc.distance_calculation(data, references, Vp, Vs)
            nevs,mevs=num.shape(distances)
            locations=references
            for evid in range(nref,nevs):
                locations=molquake_loc.dgs_locator(evid, locations, distances)
            multiloc.append(locations)
        multiloc=num.array(multiloc)
        return multiloc

    def uncertainty_estimation_2sta(self, data, references, niters, Vp, Vps=num.sqrt(3)): #not yet complete
        nref,mref=num.shape(references)
        multiloc=[]
        if len(Vp)!=2:
            print('Give a list with only 2 elements for Vp!')
        Vpmin=num.min(Vp); Vpmax=num.max(Vp)
        for niter in range(niters): #this can be parallelized
            Vp=Vpmin+(2*num.random.random_sample()-1)*(Vpmax-Vpmin)
            Vs=Vp/Vps
            distances=molquake_loc.distance_calculation_2stations(data, references, Vp, Vs)
            nevs,mevs=num.shape(distances)
            locations=references
            for evid in range(nref,nevs):
                locations=molquake_loc.dgs_locator(evid, locations, distances)
            multiloc.append(locations)
        multiloc=num.array(multiloc)
        return multiloc

    def store_results(self,out_path,filename):
        f_out=os.path.join(out_path,filename)
        try:
          num.savez(f_out+'.molres', locations=self.locations)
        except:
          print('Wrong path!!! I will use the current directory to store the files!!')
          num.savez(filename+'.molres', locations=self.locations)

    def catalogue_creation(self, evids, origin, nref, filename, output_path):
        fout=os.path.join(output_path,filename)
        nev,mev=num.shape(self.locations)
        print('number of located events',nev)
        with open(filename,'w') as f:
            f.write('Id Lat Lon Depth \n')
            for i in range(nev):
                lat,lon=LatLongUTMconversion.UTMtoLL(23, self.locations[i,1]+origin[1], self.locations[i,0]+origin[0],origin[2])
                if evids[i]!='REFERENCE':
                    event=evids[i].strftime('%Y%m%dT%H%M%S.%f')
                    f.write(event +' '+str(lat)+' '+str(lon)+' '+str(self.locations[i,2]/1000)+'\n')
