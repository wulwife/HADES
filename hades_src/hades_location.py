import numpy as num
import datetime
import LatLongUTMconversion
import os
import sys
import matplotlib.pyplot as plt

km=1000.

class hades_location:

    def __init__(self, input_obj, output_path):
        self.input=input_obj
        self.output_path=output_path

    def location(self, filename):
        distances=(self.input).distances
        references=(self.input).references
        nref,mref=num.shape(references)
        nevs,mevs=num.shape(distances)
        for i_ev in range(nref,nevs):
            sys.stdout.write(' Locating events %3d %% \r' %((i_ev/nevs)*100))
            references=hades_location.__dgslocator(i_ev, references, distances)
            sys.stdout.flush()
        self.locations=references
        self.__catalogue_creation(filename)
        self.__plot_results(filename)
        sys.stdout.write('\n')

    def __dgslocator(event, references, distances):
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
                #dij=distances[i,j]
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

    def __catalogue_creation(self, filename):
        fout=os.path.join(self.output_path,filename)
        nev,mev=num.shape(self.locations)
        evids=(self.input).events
        print('Location process completed, number of located events: %d '%(nev))
        catalogue=[]
        with open(filename+'.txt','w') as f:
            f.write('Id Lat Lon Depth Station(s) Ts-Tp\n')
            for i in range(nev):
                lat,lon=LatLongUTMconversion.UTMtoLL(23, self.locations[i,1]+(self.input).origin[1], self.locations[i,0]+(self.input).origin[0],(self.input).origin[2])
                depth=self.locations[i,2]/1000
                event=evids[i]
                t_string=' '
                for sta in (self.input).sel_sta:
                    if sta in (self.input).data[event].keys():
                        tsp=(self.input).data[event][sta][-1]
                        t_string=t_string+sta+' %5.3f '%(tsp)
                f.write(event+' '+'%6.4f '%(lat)+' '+'%6.4f '%(lon)+' '+'%3.1f '%(depth)+' '+t_string+'\n')
                catalogue.append([event,lat,lon,depth])
        self.catalogue=num.array(catalogue)


    def __plot_results(self, filename):
        c1='#4285F4'
        c2='#EA4335'
        c4='#34A853'
        c3='#FBBC05'
        nref=len((self.input).refevid)
        fig=plt.figure(figsize=(10.0,15.0))
        ax1=plt.subplot(111)
        ax1.scatter(self.locations[nref:,0],self.locations[nref:,1], s=50, c=c1)
        ax1.scatter(self.locations[0:nref,0],self.locations[0:nref:,1],s=50, c=c3)
        for sta in (self.input).stations.keys():
            station=(self.input).stations[sta]
            ax1.scatter(station[0], station[1], c=c2, marker='v', s=200, zorder=3, linewidth=0.5)
        num.min
        ax1.grid('on')
        ax1.set_aspect('equal')
        fout=os.path.join(self.output_path,filename)
        plt.savefig(fout+'.eps')
