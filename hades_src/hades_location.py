import numpy as num
import datetime
import latlon2cart
import os
import sys
import matplotlib.pyplot as plt

km=1000.

class hades_location:


    def __init__(self, input_obj, output_path):
        self.input=input_obj
        self.output_path=output_path
        
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    def location(self, filename, joint=False, fixed_dist=True, anchor=None):
        # joint              : bool, replace the incremental dgslocator loop
        #                     with a single classical-MDS solve on the full
        #                     distance matrix (masters included). Requires
        #                     distance_calculation(..., refmode='tsp') or the
        #                     master-master distances are still frozen to the
        #                     catalogue and this has no effect on them.
        #                     the referenece is the first index element in
        #                     references
        # anchor              : int or None. If given, after the above, the
        #                     whole located cluster is rigidly rotated (and
        #                     optionally mirrored) about a SINGLE reference
        #                     event (index `anchor` among the masters) to
        #                     best fit the observed ts-tp at every station
        #                     simultaneously. See __anchor_and_orient below.
        # None of this runs unless explicitly requested: with the defaults
        # (update_references=0, joint=False, anchor=None) this method is
        # byte-for-byte the original location().
        references=(self.input).references
        distances=(self.input).distances
        nref,mref=num.shape(references)
        nevs,mevs=num.shape(distances)
        nmaster=num.shape(references)[0]
        if joint:                                                       
            self.locations=self.__joint_relocation(nmaster, references, distances)
            self.__anchor_and_orient(0, demean=True, mirror=True, flip_z=False)
        else:
            for i_ev in range(nref,nevs):
                sys.stdout.write(' Locating events %3d %% \r' %((i_ev/nevs)*100))
                references=hades_location.__dgslocator(i_ev, references, distances, fixed_dist)
                sys.stdout.flush()
            self.locations=references

        if anchor is not None:
            self.__anchor_and_orient(anchor, demean=True, mirror=True, flip_z=False)
        self.__catalogue_creation(filename)
        self.__catalogue_creation_cartesian(filename)
        self.__plot_results(filename)
        sys.stdout.write('\n')


    def __dgslocator(event, references, distances, fixed_dist=True):
        '''event is the id of the event you want locate
        references is an object array of the form ['eventid',x,y,z]
        this method returns the event location and the updated the reference locations
        that include the new event located
        fixed distances means that the inter-event distances will not change during the relocation
        process (i.e. the input distance matrix remain freezed). if fixed_dist is false at each reloaction
        new inter-event distances will be recalculated from the new reloacted events'''

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
                if fixed_dist:
                    dij=distances[i,j]
                else:
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

    def __joint_relocation(self, nmaster, references, distances):
        '''Locates ALL the events - masters included - with a single
        classical (Torgerson) MDS on the full n x n distance matrix, instead
        of the incremental dgslocator loop.

        Classical scaling honours every one of the n(n-1)/2 pairwise
        distances at once (it minimises the Frobenius norm between the
        observed and the reproduced Gram matrix), so this is the
        maximum-internal-consistency solution for the whole dataset: no
        event, master included, is treated differently from any other.

        The MDS configuration is only defined up to an arbitrary
        translation, rotation and reflection. That gauge freedom is removed
        by rigidly fitting (Procrustes, no rescaling) the nmaster embedded
        master points onto the master reference coordinates (`references`,
        i.e. self.input.references or self.input.rel_references), and
        applying the SAME rigid transform to every event. The residual left
        on each master after that fit is not solver noise: it is the part
        of the mismatch between the a-priori master coordinates and the
        observed ts-tp pattern that no rigid motion can remove.

        Requires distance_calculation(..., refmode='tsp'): otherwise the
        master-master block of `distances` is still the catalogue distances
        and the masters will simply be reproduced exactly, hiding any real
        inconsistency instead of showing it.
        '''
        D=distances
        n=num.shape(D)[0]

        J=num.eye(n)-num.ones((n,n))/n
        B=-0.5*num.dot(J,num.dot(D**2,J))
        B=(B+B.T)/2.
        w,v=num.linalg.eigh(B)
        order=num.argsort(w)[::-1]
        w=w[order]; v=v[:,order]
        neg=num.sum(w[0:3]<0)
        w3=num.clip(w[0:3],0,None)
        if num.sum(w3>0)<3 or (w[2]>0 and w[0]/max(w[2],1e-12)>1e4):
            sys.stdout.write(' [warn] distance matrix close to rank 2: with two '
                             'stations only two directions per event are '
                             'independently observed, the third comes only from '
                             'the master alignment below\n')
        if neg>0:
            sys.stdout.write(' [warn] %d negative eigenvalue(s) in the MDS solution '
                             '(non-Euclidean distance matrix - noisy/inconsistent '
                             'picks); clipped to zero\n'%neg)
        Y=v[:,0:3]*num.sqrt(w3)

        X=references[0:nmaster,:]
        Ym=Y[0:nmaster,:]
        XC=num.mean(X,axis=0); YC=num.mean(Ym,axis=0)
        U1,S1,V1=num.linalg.svd(num.dot((X-XC).T,(Ym-YC)),full_matrices=True)
        Q=num.dot(U1,V1)
        Yall=num.dot(Q,(Y-YC).T).T+XC

        resid=num.sqrt(num.sum((Yall[0:nmaster,:]-X)**2,axis=1))
        for k in range(nmaster):
            sys.stdout.write(' Master %d residual after joint fit: %8.1f m\n'%(k,resid[k]))
        sys.stdout.write(' Joint relocation: %d events, master residual rms %8.1f m, max %8.1f m\n'
                         %(n,num.sqrt(num.mean(resid**2)),num.max(resid)))
        sys.stdout.flush()
        return Yall


    def __anchor_and_orient(self, anchor, demean=True, mirror=True, flip_z=False):
        '''Fixes the absolute position of the cluster on a SINGLE reference
        event (`anchor`, an index among the masters) and determines the
        remaining degrees of freedom - azimuthal rotation about that event,
        and the mirror ambiguity of the two-station geometry - by fitting the
        observed ts-tp at ALL the selected stations simultaneously.

        Meant to run after update_references or joint, where the masters are
        no longer individually pinned to absolute coordinates: at that point
        the only quantities still tied to absolute space are the position of
        the anchor and the orientation that reproduces the observations.

        demean    fit the pattern of ts-tp instead of its absolute value
                  (recommended unless kv and the station elevations are
                  trusted in an absolute sense).
        mirror    also test the reflection of the horizontal plane.
        flip_z    also test the reflection of the vertical axis.
        '''
        ntheta=721
        Vp=(self.input).vp; Vs=(self.input).vs
        kv=(Vp*Vs)/(Vp-Vs)
        depthref=(self.input).origin[-1]
        stalist=[sta for sta in (self.input).sel_sta if sta in (self.input).stations]
        if len(stalist)==0:
            sys.stdout.write(' Anchoring skipped: no station coordinates available\n')
            return

        evids=(self.input).events
        obs={}; msk={}
        for sta in stalist:
            o=num.array([(self.input).data[ev][sta][-1]
                         if sta in (self.input).data[ev] else num.nan for ev in evids])
            obs[sta]=o
            msk[sta]=num.isfinite(o)

        X=num.copy(self.locations)
        X[:,2]=X[:,2]+depthref
        target=num.array((self.input).references[anchor], dtype=float)

        def transform(my,mz,theta):
            c=((X[:,0]-X[anchor,0])+1j*my*(X[:,1]-X[anchor,1]))*num.exp(-1j*theta)
            Y=num.zeros(num.shape(X))
            Y[:,0]=c.real+target[0]
            Y[:,1]=c.imag+target[1]
            Y[:,2]=mz*(X[:,2]-X[anchor,2])+target[2]
            return Y

        def misfit(Y):
            num_=0.; den=0
            for sta in stalist:
                m=msk[sta]
                s=(self.input).stations[sta]
                d=num.sqrt((Y[m,0]-s[0])**2+(Y[m,1]-s[1])**2+(Y[m,2]-s[2])**2)
                calc=d/kv; o=obs[sta][m]
                if demean:
                    calc=calc-num.mean(calc); o=o-num.mean(o)
                num_+=num.sum((calc-o)**2); den+=num.size(o)
            return num.sqrt(num_/den)

        thetas=num.linspace(0,2*num.pi,ntheta)
        best=(None,None,None,num.inf)
        for my in ([1,-1] if mirror else [1]):
            for mz in ([1,-1] if flip_z else [1]):
                rms=num.array([misfit(transform(my,mz,t)) for t in thetas])
                k=int(num.argmin(rms))
                lo=thetas[max(k-1,0)]; hi=thetas[min(k+1,ntheta-1)]
                fine=num.linspace(lo,hi,101)
                rmsf=num.array([misfit(transform(my,mz,t)) for t in fine])
                kf=int(num.argmin(rmsf))
                sys.stdout.write(' Anchor fit: mirror_y %+d mirror_z %+d -> theta %6.1f deg, '
                                 'rms %8.4f s\n'%(my,mz,num.degrees(fine[kf]),rmsf[kf]))
                if rmsf[kf]<best[3]:
                    best=(my,mz,fine[kf],rmsf[kf])
        my,mz,theta,rms=best
        self.locations=transform(my,mz,theta)
        self.locations[:,2]=self.locations[:,2]-depthref
        self.anchor_solution={'anchor':(self.input).events[anchor],'mirror_y':my,
                              'mirror_z':mz,'theta_deg':num.degrees(theta),'rms':rms}
        sys.stdout.write(' Cluster anchored on %s: mirror_y %+d, mirror_z %+d, '
                         'theta %.1f deg, rms %.4f s\n'
                         %((self.input).events[anchor],my,mz,num.degrees(theta),rms))
        sys.stdout.flush()

    def __absolute_cluster_location(self,filename):
        #currently only search along strike is implemented
        Vp=(self.input).vp
        Vs=(self.input).vs
        kv=(Vp*Vs)/(Vp-Vs)
        stations=(self.input).stations
        depth=(self.input).origin[-1]
        thetas=num.arange(0,41)*0.025*num.pi*2
        evtsps=self.__initialize_tsp_db(stations)
        rms_min=1E10
        zrot=self.locations[:,2]+depth
        for ysign in [-1,1]:
            for theta in thetas:
                crot=(self.locations[:,0]+1j*ysign*self.locations[:,1])*num.exp(-1j*theta)
                rms=self.__rms_theta_calculation(crot.real,crot.imag,zrot,evtsps,kv,stations)
                if rms < rms_min:
                    rms_min=rms
                    theta_best=theta
                    ysign_best=ysign
        crot=(self.locations[:,0]+1j*ysign_best*self.locations[:,1])*num.exp(-1j*theta_best)
        self.locations[:,0]=crot.real
        self.locations[:,1]=crot.imag
        self.locations[:,2]=zrot
        pca_max=0
        for theta in thetas:
            crot=(self.locations[:,0]+1j*ysign_best*self.locations[:,1])*num.exp(-1j*theta)
            pca=self.__pca_theta_calculation(crot.real,crot.imag,zrot,evtsps,stations)
            if pca > pca_max:
                pca_max=pca
                theta_best=theta
        crot=(self.locations[:,0]+1j*self.locations[:,1])*num.exp(-1j*theta_best)
        self.locations[:,0]=crot.real
        self.locations[:,1]=crot.imag
        self.locations[:,2]=zrot
        self.__catalogue_creation(filename)
        self.__plot_results(filename)


    def __rms_theta_calculation(self,xobs,yobs,zobs,evtsps,kv,stations):
        rms=0
        for sta in stations.keys():
            dx=(xobs-(self.input).stations[sta][0])
            dy=(yobs-(self.input).stations[sta][1])
            dz=(zobs-(self.input).stations[sta][2])
            tsp_obs=num.array(evtsps[sta])
            tsp_obs=tsp_obs-num.mean(tsp_obs)
            tsp_calc=num.sqrt(dx**2+dy**2+dz**2)/kv
            tsp_calc=tsp_calc-num.mean(tsp_calc)
            rms+=num.sqrt(num.sum((tsp_calc-tsp_obs)**2)/num.size(tsp_obs))
        rms=rms/len(stations.keys())
        return rms

    def __pca_theta_calculation(self,xobs,yobs,zobs,evtsps,stations):
        rect=1
        signs=[]
        for sta in stations.keys():
            X=num.zeros([num.size(xobs),2])
            dx=(xobs-(self.input).stations[sta][0])
            dy=(yobs-(self.input).stations[sta][1])
            dz=(zobs-(self.input).stations[sta][2])
            tsp=num.array(evtsps[sta])
            #tsp=tsp_obs-num.mean(tsp_obs)
            dist=num.sqrt(dx**2+dy**2+dz**2)
            ir_dist=num.argsort(dist)
            X[:,0]=dist[ir_dist]
            X[:,1]=tsp[ir_dist]
            M=num.mean(X.T, axis=1)
            C=X-M
            V=num.cov(C.T)
            values, vectors = num.linalg.eigh(V)
            sign=num.sign(vectors[0,1]*vectors[1,1])
            signs.append(1*sign)
            rect=rect*(num.max(values)/num.min(values))
        if signs[0]>0 and signs[1]>0:
            rect=rect
        else:
            rect=-1*rect
        return rect

    def __initialize_tsp_db(self,stations):
        evtsps={}
        for sta in stations.keys():
            evtsps[sta]=[]
            for event in (self.input).events:
                evtsps[sta].append((self.input).data[event][sta][-1])
        return evtsps

    def __cluster_orientation(self):
        #currently only search along strike is implemented
        ref=(self.input).references
        ref1_x=ref[1,0]-ref[0,0]
        ref1_y=ref[1,1]-ref[0,1]
        theta1=num.atan2(ref1_y,ref1_x)
        #theta2=num.atan2(ref2_y,ref2_x)

    def __catalogue_creation(self, filename):
        fout=os.path.join(self.output_path,filename)
        nev,mev=num.shape(self.locations)
        evids=(self.input).events
        latref,lonref=(self.input).origin[0],(self.input).origin[1]
        orig=latlon2cart.Coordinates(latref,lonref,0)
        print('Location process completed, number of located events: %d '%(nev))
        catalogue=[]
        with open(fout+'.txt','w') as f:
            f.write('Id Lat Lon Depth Station(s) Tp Ts-Tp\n')
            for i in range(nev):
                lat,lon,_=orig.cart2geo(self.locations[i,0],self.locations[i,1],0)
                depth=self.locations[i,2]/1000
                event=evids[i]
                t_string=' '
                for sta in (self.input).sel_sta:
                    if sta in (self.input).data[event].keys():
                        tsp=(self.input).data[event][sta][-1]
                        tid=(self.input).data[event][sta][0]
                        t_string=t_string+sta+' %5.3f '%(tsp)+str(tid)+ ' '
                f.write(event+' '+'%6.4f '%(lat)+' '+'%6.4f '%(lon)+' '+'%3.1f '%(depth)+' '+t_string+'\n')
                catalogue.append([event,lat,lon,depth])
        self.catalogue=num.array(catalogue)

    def __catalogue_creation_cartesian(self, filename):
        fout=os.path.join(self.output_path,filename)
        nev,mev=num.shape(self.locations)
        evids=(self.input).events
        latref,lonref=(self.input).origin[0],(self.input).origin[1]
        print('Catalogues generated')
        catalogue_cart=[]
        with open(fout+'cartesian.txt','w') as f:
            f.write('Id X Y Z Station(s) Tp Ts-Tp\n')
            f.write('Ref Lat :' + str(latref) + ', Ref Lon :'+ str(lonref) + '\n')
            for i in range(nev):
                x,y=self.locations[i,0]/1000,self.locations[i,1]/1000
                depth=self.locations[i,2]/1000
                event=evids[i]
                t_string=' '
                for sta in (self.input).sel_sta:
                    if sta in (self.input).data[event].keys():
                        tsp=(self.input).data[event][sta][-1]
                        tid=(self.input).data[event][sta][0]
                        t_string=t_string+sta+' %5.3f '%(tsp)+str(tid)+ ' '
                f.write(event+' '+'%6.3f '%(x)+' '+'%6.3f '%(y)+' '+'%3.1f '%(depth)+' '+t_string+'\n')
                catalogue_cart.append([event,x,y,depth])
        self.catalogue_cart=num.array(catalogue_cart)


    def __plot_results(self, filename):
        c1 = '#4285F4'
        c2 = '#EA4335'
        c4 = '#34A853'
        c3 = '#FBBC05'
        nref = len((self.input).refevid)

        # event coordinates
        ex = self.locations[:, 0]
        ey = self.locations[:, 1]
        ez = self.locations[:, 2]

        # station coordinates
        stx, sty, stz = [], [], []
        for sta in (self.input).stations.keys():
            station = (self.input).stations[sta]
            stx.append(station[0])
            sty.append(station[1])
            stz.append(station[2])
        stx = num.array(stx)
        sty = num.array(sty)
        stz = num.array(stz)

        # common span for all axes -> identical, comparable panels
        pad = 0.05
        lims = []
        for v in (ex,ey,ez):
            lims.append((v.min(), v.max()))
        span = max([b - a for a, b in lims])
        span = span * (1.0 + 2.0 * pad) if span > 0 else 1.0
        xlim, ylim, zlim = [(0.5 * (a + b) - 0.5 * span,
                             0.5 * (a + b) + 0.5 * span) for a, b in lims]

        fig = plt.figure(figsize=(12.0, 12.0))

        # --- map view: X-Y ---
        ax1 = plt.subplot(221)
        ax1.scatter(ex[nref:], ey[nref:], s=50, c=c1)
        ax1.scatter(ex[0:nref], ey[0:nref], s=50, c=c3)
        ax1.scatter(stx, sty, c=c2, marker='v', s=200, zorder=3, linewidth=0.5)
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')

        # --- cross-section: Z-Y (view from the side, depth on x) ---
        ax2 = plt.subplot(222, sharey=ax1)
        ax2.scatter(ez[nref:], ey[nref:], s=50, c=c1)
        ax2.scatter(ez[0:nref], ey[0:nref], s=50, c=c3)
        ax2.set_xlim(zlim)
        ax2.set_ylim(ylim)
        ax2.set_xlabel('Z')
        ax2.set_ylabel('Y')

        # --- cross-section: X-Z (depth on y) ---
        ax3 = plt.subplot(223, sharex=ax1)
        ax3.scatter(ex[nref:], ez[nref:], s=50, c=c1)
        ax3.scatter(ex[0:nref], ez[0:nref], s=50, c=c3)
        ax3.set_xlim(xlim)
        ax3.set_ylim(zlim)
        ax3.set_xlabel('X')
        ax3.set_ylabel('Z')

        for ax in (ax1, ax2, ax3):
            ax.grid('on')
            ax.set_aspect('equal', adjustable='box')

        plt.subplot(224).axis('off')
        fig.tight_layout()

        fout = os.path.join(self.output_path, filename)
        plt.savefig(fout + '.pdf')
        plt.close(fig)