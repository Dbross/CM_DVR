#!/usr/bin/env python3
""" This program was written following 
Citation: The Journal of Chemical Physics 96, 1982 (1992); doi: 10.1063/1.462100
as a guide
Can use  scipy.sparse.linalg.eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')[source] as a way of doing Lanczos...."""
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (bytes, str, open, super, range, zip, round, input, int, pow, object)

def constants(CODATA_year=2010):
    """ CODATA constants used in program and other global defintions including number of points and numpy datatype for reals"""
    global numpy_precision, num_points, num_print_digits
    numpy_precision="np.float64"
    num_points=101
    num_print_digits=3
    global planckconstant, light_speed, Rydberg, electron_charge, amu, bohr, e_mass, hartreetocm
    light_speed= 299792458 # m/s
    if CODATA_year==2010:
        planckconstant=6.62606957E-34 # Js
        amu=1.660538921E-27 #kg
        e_mass= 9.10938215E-31 #kg
        electron_charge=1.602176565E-19  #C
        rydberg= 109737.31568539 # cm-1
    elif CODATA_year==2006:
        planckconstant= 6.626068963E-34 # Js
        amu= 1.660538782E-27 # kg
        e_mass= 9.10938291E-31 #kg
        electron_charge= 1.602176487E-19 # C
        rydberg= 109737.31568527 #cm-1
    else:
        from sys import exit
        exit('Constants not found')
    bohr= (float(5) * light_speed * electron_charge**2)/ (planckconstant*rydberg) 
    hartreetocm= 2*rydberg
#    hartreetocm=219474.6313717

class potential:
    """ Define this class, along with keywords"""
# anyline starting with ! or # is commented out
    def __init__(self):
        self.r=[]
        self.energy=[]
        self.coordtypes=[]
        self.mass=[]
        self.pts=[]
        self.mingrid=False
        self.harmonicfreq=[]
        self.minpos=[]
        self.fiteignum=[]
        self.fiteigval=[]
        self.massadjust=False
        self.saveeigrequested=False
        self.preplot=False
        self.plotit=False
        self.print2ndderiv=False
        self.printeigenval=True
        self.useprevpoint=False
        self.prevpoints=[]
        self.petsc=False
        self.printpetsc=False

    def readpotential(self,inp):
        #potential should have coordinate and units as main input
        commentoutchars=['!','#']
        types=['angular_2pi','radial','angular_pi']
        lines=openandread(inp)
        r_unitconversion=1.0
        r=[]
#       r_unitconversion=1.88972613
        import re
        numeric_const_pattern = r"""[-+]?(?: (?: \d* \. \d+ ) | (?: \d+ \.? ) ) (?: [Ee] [+-]? \d+ ) ?"""
        rx=re.compile(numeric_const_pattern,re.VERBOSE)
        self.emin=0
        massassigned=False
        for x in lines:
            if x[0] not in commentoutchars:
                if '---' in x.lower():
                    break
                elif 'mass' in x.lower() and not massassigned:
                    self.mass=rx.findall(x)
                    for x in range(len(self.mass)):
                        self.mass[x]=float(self.mass[x])
                    print('Input reduced mass of {0} amu.'.format(self.mass))
                    massassigned=True
                elif 'mass' in x.lower() and massassigned:
                    from sys import exit
                    exit('mass defined in potential twice')
                elif 'bohr' in x.lower():
                    print('reading potential as bohr, this should only be set once.')
                    r_unitconversion=1.0
                elif 'emin' in x.lower():
                    self.emin=rx.findall(x)
                    if len(self.emin)>1:
                        for x in self.emin:
                            x=float(x)
                        print('found multiple energy minimum {0} using {1} as minimum'.format(self.emin,min(self.emin)))
                        self.emin=min(self.emin)
                    else:
                        self.emin=float(self.emin[0])
                    print('shifting energy minimum to {} a.u.'.format(self.emin))
                elif types[0] in x.lower() or types[1] in x.lower() or types[2] in x.lower():
                    typelist=x.lower().replace(',',' ').split()
                    for y in typelist:
                        if y in types:
                            self.coordtypes.append(y)
                elif 'preplot' in x.lower():
                    self.preplot=True
                    self.preplotval=re.findall(r'\b\d+\b', x)
                    if len(self.preplotval)>=1:
                        for x in range(len(self.preplotval)):
                            self.preplotval[x]=int(self.preplotval[x])
                    else:
                        self.preplotval=[int(2000)]
                elif 'plotit' in x.lower():
                    self.plotit=True
                elif 'petsc' in x.lower():
                    self.petsc=True
                    self.printpetsc=True
                    self.numsol=re.findall(r'\b\d+\b', x)
                    if len(self.numsol)>=1:
                        self.numsol=int(self.numsol[0])
                    else:
                        self.numsol=6
                elif 'mingrid' in x.lower():
                    self.mingrid=True
                elif 'printderiv' in x.lower():
                    self.print2ndderiv=True
                elif 'noprinteigval' in x.lower():
                    self.printeigenval=False
                elif 'fiteigval' in x.lower():
                    self.fiteigval=rx.findall(x)
                    for x in range(len(self.fiteigval)):
                        self.fiteigval[x]=float(self.fiteigval[x])/ hartreetocm
                elif 'fiteignum' in x.lower():
                    self.fiteignum=re.findall(r'\b\d+\b', x)
                    for x in range(len(self.fiteignum)):
                        self.fiteignum[x]=int(self.fiteignum[x])
                    self.petsc=True
                    self.numsol=max(self.fiteignum)+1
                elif 'minpos' in x.lower():
                    self.minpos=rx.findall(x)
                    print('Will use {0} as starting minimum for 2nd derivative evaluation if requested.'.format(self.minpos))
                    for x in range(len(self.minpos)):
                        self.minpos[x]=float(self.minpos[x])
                elif 'harmonic' in x.lower():
                    self.harmonicfreq=rx.findall(x)
                    print('Will adjust reduced mass to match harmonic frequencies of {0} cm-1.'.format(self.harmonicfreq))
                    for x in range(len(self.harmonicfreq)):
                        self.harmonicfreq[x]=float(self.harmonicfreq[x])/hartreetocm
                elif 'angstrom' in x.lower():
                    print('reading potential as angstrom, this should only be set once.')
                    r_unitconversion=(1.0/bohr)
                else:
                    linesplit=x.split()
                    rtmp=[]
                    for x in range(0,len(linesplit)-1):
                        rtmp.append(float(linesplit[x])*r_unitconversion)
                    r.append(rtmp)
                    self.energy.append(float(linesplit[len(linesplit)-1]))
                    if len(r[-1])!=len(self.coordtypes):
                        from sys import exit
                        print(x)
                        exit('number of coordinates given inconsistent with coordinate types given')
            else:
                print('{0} commented out'.format(x[1:].strip()))
        for i in range(len(r[0])):
            rtmp=[r[j][i] for j in range(len(r))]
            self.r.append(rtmp)
        if len(self.fiteignum)!=len(self.fiteigval):
            print('{0} eigenvalues to be fit given and {1} assignments given for them'.format(len(self.fiteigval),len(self.fiteignum)))
            from sys import exit
            exit()
        if len(self.mass)!=len(self.coordtypes):
            print('{0} masses given and {1} coordinate types give'.format(len(self.mass),len(self.coordtypes)))
            from sys import exit
            exit()
        import numpy as np
        self.r=np.array(self.r,dtype=eval(numpy_precision))
        self.energy=np.subtract(np.array(self.energy,dtype=eval(numpy_precision)),self.emin)
        """ Radial terminates with PiB walls; angular_2pi repeats; angular_pi terminates with PiB walls at 0 and pi"""
        coordtypedict={'radial': 'r', 'angular_2pi': 'phi', 'angular_pi': 'theta'}
        for x in range(len(self.coordtypes)):
            self.coordtypes[x]=coordtypedict[self.coordtypes[x]]
        for x in range(len(self.r)):
            self.pts.append(len(np.unique(self.r[x])))

    def xlsx(self):
        import xlsxwriter
        workbook = xlsxwriter.Workbook('tmp.xlsx')
        worksheet = workbook.add_worksheet()
        for i in range(len(self.r)):
            for j in range(len(self.r[i])):
                worksheet.write(j,i,self.r[i][j])
                if i==0:
                    worksheet.write(j,len(self.r),self.energy[j])
        workbook.close()

    def solve(self):
        if len(self.coordtypes)==1:
#            from timeit import Timer
#            t = Timer(lambda: spline1dpot(pts,mass,coordtypes,Energies,r))
#            print('time={0}'.format(t.timeit(number=10)))
            self.spline1dpot(self.pts,self.mass,self.coordtypes,self.energy,self.r)
        elif len(self.coordtypes)==2:
#            from timeit import Timer
#            t = Timer(lambda: spline2dpot(pts,mass,coordtypes,Energies,r))
#            print('time={0}'.format(t.timeit(number=1)))
            import numpy as np
            self.spline2dpot([ x for x in self.pts ],[ x for x in self.mass ],self.coordtypes,self.energy,self.r)
        else:
            raise ValueError("not implemented for %d dimensions" % (len(self.coordtypes)))
            Ham=H_array(pts=pts,mass=mass,V=Energies,qmax=np.amax(r,axis=1),qmin=np.amin(r,axis=1),coordtype=coordtypes)
            eigenval, eigenvec=np.linalg.eigh(Ham)
            Esort=np.sort(eigenval*hartreetocm)
            Etoprint=int(len(Esort)/2)
            for x in range(Etoprint):
                print('{0:.{1}f}'.format(round(Esort[x],num_print_digits),num_print_digits))

    def fitfundamental(self):
        from scipy.optimize import minimize
        import numpy as np
        c=[]
        for x in range(len(self.fiteignum)):
            c.append(self.mass[x])
#            c.append(np.multiply(self.mass[x],\
#                    np.square(np.divide(np.subtract(self.eigenval[self.fiteignum[x]],self.eigenval[0]),self.fiteigval[x]))))
        result=minimize(self.calcfreqminusactualfreq,c,method='L-BFGS-B',options={'ftol':1e-04,'gtol':1e-04,'eps':1e-07})
        print(result)
        print('Final Eigenvalues')
        for x in range(len(self.fiteignum)): 
            print(self.eigenval[self.fiteignum],np.subtract(self.eigenval[self.fiteignum[x]],self.eigenval[0]))
        print('Final Mass {0}'.format(self.mass))

    def calcfreqminusactualfreq(self,c):
        self.mass=c
        print(self.mass)
        self.solve()
        tot=0.0
        from numpy import add, subtract, abs
        for x in range(len(self.fiteignum)):
            tot=add(abs(subtract(self.fiteigval[x],subtract(self.eigenval[self.fiteignum[x]],self.eigenval[0]))),tot)
        return tot
    
    def printeigenvals(self):
        from numpy import multiply
        Esort=multiply(self.eigenval,hartreetocm)
        Etoprint=int(len(Esort)/2)
        for x in range(Etoprint):
            print('{0:.{1}f}'.format(round(Esort[x],num_print_digits),num_print_digits))

    def spline1dpot(self,pts,mass,coordtypes,Energies_raw,r_raw):
        import numpy as np
        xmin,emin=return1dsplinemin(r_raw[0],Energies_raw)
        #r=r_raw-np.min(xmin)
        r=r_raw
        Energies=Energies_raw-np.min(emin)
        from scipy.interpolate import splrep
        Ener_spline=splrep(r[0],Energies,s=0)
        from scipy.interpolate import splev
        if self.print2ndderiv:
            print('{0:10} {1:15} {2:15}'.format('Position','2nd Derivative','Eval cm-1'))
            smooth_spline=splrep(r[0],Energies,s=1e-7)
            Energ=splev(xmin,smooth_spline)
            derivs=splev(xmin,smooth_spline,der=2)
            for x in range(len(xmin)):
                print('({0:10.4e} {1:15.4e} {2:15.2f})'.format(xmin[x],derivs[x],Energ[x]*hartreetocm))
        if len(self.harmonicfreq)==1 and not self.massadjust:
#            raise ValueError("not implemented for %d dimensions" % (len(self.coordtypes)))
            if len(self.minpos)==0:
                minpos=(min(r[0])+max(r[0]))/2.0
            else:
                minpos=self.minpos[0]
            mass=np.multiply(np.divide(splev(np.array(minpos),Ener_spline,der=2),np.power(self.harmonicfreq,2)),np.divide(e_mass,amu)) 
            print('Adjusted potential to use mass of {0} based on harmonic frequency.'.format(mass))
            self.massadjust=True
        xnew = np.linspace(min(r[0]),max(r[0]), num=num_points)
        vfit= splev(xnew, Ener_spline, der=0)
        Ham=H_array(pts=pts,mass=mass,V=Energies,qmax=np.amax(r,axis=1),qmin=np.amin(r,axis=1),coordtype=coordtypes)
        eigenval, eigenvec=np.linalg.eigh(Ham)
        eindex=np.argsort(eigenval)
        self.eigenval, eigenvec= eigenval[eindex], np.transpose(eigenvec[:,eindex])
        Esort=(self.eigenval*hartreetocm)
        Etoprint=int(len(Esort)/2)
        maxpot=np.max(vfit)*hartreetocm
        i=0
        for i in range(Etoprint):
            if Esort[i]>maxpot:
                Etoprint=i-1
                break
        if self.plotit:
            vfitcm=vfit*hartreetocm
            import matplotlib.pyplot as plt
            if self.print2ndderiv:
                der1= splev(xnew, smooth_spline, der=1)
                der2= splev(xnew, smooth_spline, der=2)
                plot1dfit(self.r[0],self.energy,xnew,der1,title='Spline 1st deriv',block=False)
                plot1dfit(self.r[0],self.energy,xnew,der2,title='Spline 2nd deriv',block=False)
            plt.figure()
            plt.plot(r[0],Energies*hartreetocm,linestyle='none',marker='o')
            plt.plot(xnew,vfitcm,linestyle='solid',marker='None')
            mincut=[]
            maxcut=[]
            maxpot=np.max(np.multiply(Energies,hartreetocm))
            for x in range(Etoprint):
                for y in range(len(xnew)):
                    if Esort[x]>vfitcm[y]:
                        mincut.append(xnew[y])
                        break
            for x in range(Etoprint):
                for y in range(len(xnew)-1,1,-1):
                    if Esort[x]>vfitcm[y]:
                        maxcut.append(xnew[y])
                        break
            for x in range(Etoprint):
                plt.plot((mincut[x],maxcut[x]),(Esort[x],Esort[x]),linestyle='solid')
            plt.legend(['Points', 'Cubic Spline'])
            plt.title('Cubic-spline interpolation')
            plt.axis()
            plt.show(block=True)
#            plt.figure()
#            plt.title('ground state')
#            plt.plot(r[0],Energies/np.max(Energies))
#            plt.plot(r[0],np.square(eigenvec[0]),marker='o')
#            plt.plot(r[0],eigenvec[0],marker='x')
#            plt.show(block=False)
#            plt.figure()
#            plt.title('v=1 state')
#            plt.plot(r[0],Energies/np.max(Energies))
#            plt.plot(r[0],eigenvec[1],marker='o')
#            plt.plot(r[0],np.square(eigenvec[1]),marker='x')
#            plt.show(block=True)

    def fit1dpot(self):
        import numpy as np
        from scipy.optimize import curve_fit, root, brentq
        import inspect
        #init=[np.max(self.energy)]* (len(inspect.getfullargspec(func1d).args)-2)
        init=[0.0001]*(len(inspect.getfullargspec(func1d).args)-2) 
        popt,pconv= curve_fit(func1d, self.r[0],self.energy,p0=(0,*init))
        xgrid=np.linspace(np.min(self.r[0]),np.max(self.r[0]),100)
        popt=np.ndarray.tolist(popt)
        plot1dfit(self.r[0],self.energy,xgrid,func1d(xgrid,*popt),title='analytic func',block=False)
        plot1dfit(self.r[0],self.energy,xgrid,func1dder1(xgrid,*popt),title='1nd Deriv, analytic func',block=False)
        plot1dfit(self.r[0],self.energy,xgrid,func1dder2(xgrid,*popt),title='2nd Deriv, analytic func',block=False)
        derivgrid=func1dder1(xgrid,*popt)
        asign = np.sign(derivgrid)
        signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
        signchange[0]=0
        roots=[]
        print('From fit potential')
        print('{0:10} {1:15} {2:15}'.format('Position','2nd Derivative','Eval cm-1'))
        for x in range(len(signchange)):
            if signchange[x]==1:# and np.greater(func1dder2(np.squeeze(xgrid[x]),*popt),0):
                tmproot=brentq(func1dder1,xgrid[x-1],xgrid[x],args=tuple(popt))
                print('({0:10.4e} {1:15.4e} {2:15.2f})'\
                        .format(tmproot,func1dder2(tmproot,*popt),func1d(tmproot,*popt)*hartreetocm))
#                print(derivgrid[x],derivgrid[x-1])

    def spline2dpot(self,pts,mass,coordtypes,Energies,r,saveeigen=True):
        import numpy as np
        from scipy.interpolate import griddata
        from potgen import silentmweq, roundmasstoequal
        if self.preplot:
            for x in range(len(self.preplotval)):
                if x==len(self.preplotval)-1:
                    plot2d(r[0],r[1],Energies,wavenumber=True,title='Potential Energy Contours',block=True,includegrid=False,wavenumbercutoff=self.preplotval[x])
                else:
                    plot2d(r[0],r[1],Energies,wavenumber=True,title='Potential Energy Contours',block=False,includegrid=False,wavenumbercutoff=self.preplotval[x])
        sigfigs=4
        qmin0 =np.copy(np.min(r[0]))
        qmax0 =np.copy(np.max(r[0]) )
        qmin1 =np.copy(np.min(r[1]) )
        qmax1 =np.copy(np.max(r[1]) )
        org=[qmin0,qmax0,qmin1,qmax1]
        """ rlen is (max-min) of potential, scaled to b-a by adding two more points!"""
        rlen0=np.multiply(np.subtract(qmax0,qmin0),np.divide(np.add(float(pts[0]),1.0),np.subtract(float(pts[0]),1.0)))
        rlen1=np.multiply(np.subtract(qmax1,qmin1),np.divide(np.add(float(pts[1]),1.0),np.subtract(float(pts[1]),1.0)))
        if len(self.harmonicfreq)==2 and not self.massadjust:
            from scipy.interpolate import RectBivariateSpline, bisplev, bisplrep, spleval
            cubic2dspline= RectBivariateSpline(np.unique(r[0]), np.unique(r[1]),  np.reshape(Energies,(pts[0],pts[1])),s=1e-6)
            from scipy.optimize import minimize
            minbnds=(qmin0,qmin1)
            maxbnds=(qmax0,qmax1)
            bnds=list(zip(minbnds,maxbnds))
            if len(self.minpos)==0:
                c=[ (qmin0+qmax0)/2.0, (qmin1+qmax1)/2.0]
            else:
                c=[self.minpos[0],self.minpos[1]]
            result=minimize(return2dspline,c,args=(cubic2dspline,0,0),bounds=bnds,jac=return2dsplinetotder,method='L-BFGS-B')
            x,y=result.x
            hessx=cubic2dspline.ev(x,y,dx=2) # calculate the second partial derivitive for dq0 at minimia
            hessy=cubic2dspline.ev(x,y,dy=2) # calculate the second partial derivitive for dq1 at minimia
            if (float(hessx) < 0) or (float(hessy) < 0):
                print('Issues in hessian calculation, increase smoothing value in RectBivariateSpline in potential.spline2dpot to help')
                self.plotit=True
            print('Hessians calculated at ({2:.4e}, {3:.4e}) Eval={4:.0f} cm-1 as [{0:.9e} dq0^2, {1:.9e} dq1^2].'\
                    .format(float(hessx),float(hessy),float(x),float(y),float(result.fun)*hartreetocm))
            mass[0]=np.multiply(np.divide(hessx,np.power(self.harmonicfreq[0],2)),np.divide(e_mass,amu)) 
            mass[1]=np.multiply(np.divide(hessy,np.power(self.harmonicfreq[1],2)),np.divide(e_mass,amu)) 
            print('Adjusted potential to use mass of {0} based on harmonic frequencies.'.format(mass))
            self.massadjust=True
            if self.plotit:
                x,y= np.meshgrid(np.unique(r[0]), np.unique(r[1]))
                z=return2dspline((x,y),cubic2dspline,0,0)
                plot2d(r[0],r[1],Energies,wavenumber=True,title='Potential Energy Contours',block=False)
                plot2dgrid(x,y,z,wavenumber=True,title='PES Fit',block=True)
# this plots the 2d grid, incase you'd like to see which point corresponds to which coordinate
#        import matplotlib.pyplot as plt
#        plt.figure()
#        plt.plot(r[0]*180/np.pi,r[1]*180/np.pi,linestyle='none',marker='o')
#        labels=['{0}'.format(i) for i in range(len(r[0]))]
#        for label, x,y in zip(labels,r[0,:],r[1,:]): 
#            plt.annotate(label,
#                    xy= (x,y), xytext = (-5,5),
#                    textcoords= 'offset points', ha='right', va='bottom')
##                    ,bbox=dict(boxstyle='round,pad=0.5',fc='yellow',alpha=0.5))
##                    ,arrowprops=dict(arrowstyle ='->',connectionstyle='arc3,rad=0'))
#        plt.show()
        mw1=mwspace(coordtype=coordtypes[0],rlen=rlen0,mass=mass[0],pts=pts[0])
        mw2=mwspace(coordtype=coordtypes[1],rlen=rlen1,mass=mass[1],pts=pts[1])
        if np.abs(np.subtract(mw1,mw2))>1.0E-07 or self.mingrid:
            if self.useprevpoint:
                pts[0]=self.prevpoints[0]
                pts[1]=self.prevpoints[1]
        mw1=mwspace(coordtype=coordtypes[0],rlen=rlen0,mass=mass[0],pts=pts[0])
        mw2=mwspace(coordtype=coordtypes[1],rlen=rlen1,mass=mass[1],pts=pts[1])
        if np.abs(np.subtract(mw1,mw2))>1.0E-07 or self.mingrid:
            """ Adjust potential to have equal massweighted spacing"""
            print('Mass weighting unequal, adjusting grid\n OLD: {0:.4f}-{1:.4f} pts {2} {3:.4f}-{4:.4f} pts {5}'\
                    .format(float(qmin0),float(qmax0),int(pts[0]),float(qmin1),float(qmax1),int(pts[1])))
            a=silentmweq([ [qmax0,qmin0,coordtypes[0],pts[0],mass[0]], [qmax1,qmin1,coordtypes[1],pts[1],mass[1]] ],mingrid=self.mingrid,uselowest=self.useprevpoint)
            """ It is possible to override points and grid to be fit to here... Thinking about adding a manual option but unsure why I'd do that..."""
            qmax0,qmin0,pts[0]=np.max(a[0].grid),np.min(a[0].grid),a[0].numpoints
            qmax1,qmin1,pts[1]=np.max(a[1].grid),np.min(a[1].grid),a[1].numpoints
            from sys import exit
            if np.subtract(qmax1,org[3])>0.0:
                if np.subtract(qmax1,org[3])<1E-11:
                    qmax1=org[3]
                else:
                    exit('max of coordinate 1 exceeds potential')
            if np.subtract(qmax0,org[1])>0.0:
                if np.subtract(qmax0,org[1])<1E-11:
                    qmax0=org[1]
                else:
                    exit('max of coordinate 0 exceeds potential')
            if np.subtract(org[0],qmin0)>0.0:
                if np.subtract(org[0],qmin0)<1E-11: 
                    qmin0=org[0]
                else:
                    exit('min of coordinate 0 exceeds potential')
            if np.subtract(org[2],qmin1)>0.0:
                if np.subtract(org[2],qmin1)<1E-11: 
                    qmin1=org[2]
                else:
                    exit('min of coordinate 1 exceeds potential')
            grid_x, grid_y = np.mgrid[qmin0:qmax0:pts[0]*1j,qmin1:qmax1:pts[1]*1j]
            print('NEW: {0:.4f}-{1:.4f} pts {2} {3:.4f}-{4:.4f} pts {5}'\
                    .format(qmin0,qmax0,pts[0],qmin1,qmax1,pts[1]))
            vfit= griddata(np.transpose(np.array(r)),Energies,(grid_x,grid_y),method='cubic')
            if np.any(np.isnan(vfit)):
                exit('fit potential outside bounds')
#            mass= roundmasstoequal(mass=mass,sigfigs=sigfigs,dq1=np.divide(mw1,mass[0]),dq2=np.divide(mw2,mass[1]))
#            print('using {2} sig figs of reduced mass of [{0:.{3}e}, {1:.{3}e}] amu.'.format(float(mass[0]),float(mass[1]),sigfigs,sigfigs-1))
            self.useprevpoint=True
            self.prevpoints=pts
            if self.petsc:
                eigenval=H_array_petsc(pts=pts,mass=mass,V=np.ndarray.flatten(vfit),\
                        qmax=[qmax0,qmax1],qmin=[qmin0,qmin1],coordtype=coordtypes,numeig=self.numsol\
                        ,printpetsc=self.printpetsc)
            else:
                Ham=H_array(pts=pts,mass=mass,V=np.ndarray.flatten(vfit),\
                        qmax=[qmax0,qmax1],qmin=[qmin0,qmin1],coordtype=coordtypes)
                eigenval, eigenvec=np.linalg.eigh(Ham)
        else:   
#            mass= roundmasstoequal(mass=mass,sigfigs=sigfigs,dq1=np.divide(mw1,mass[0]),dq2=np.divide(mw2,mass[1]))
#            print('using {2} sig figs of reduced mass of [{0:.{3}e}, {1:.{3}e}] amu.'.format(float(mass[0]),float(mass[1]),sigfigs,sigfigs-1))
            if self.petsc:
                eigenval=H_array_petsc(pts=pts,mass=mass,V=Energies,qmax=np.amax(r,axis=1),qmin=np.amin(r,axis=1),\
                        coordtype=coordtypes,numeig=self.numsol,printpetsc=self.printpetsc)
            else:
                Ham=H_array(pts=pts,mass=mass,V=Energies,qmax=np.amax(r,axis=1),qmin=np.amin(r,axis=1),coordtype=coordtypes)
                eigenval, eigenvec=np.linalg.eigh(Ham)
        self.eigenval=eigenval.real.astype(eval(numpy_precision))
        eindex=np.argsort(eigenval)
        if not self.petsc:
            self.eigenval, eigenvec= eigenval[eindex], np.transpose(eigenvec[:,eindex])
#        from scipy.sparse.linalg import eigs
        """ Sparse solver doesn't give speedup for computing eigenvalues of all solutions.... """
#       eigenval, eigenvec=eigs(Ham,k=int((pts[0]*pts[1])-2),sigma=0,M=None,which='LM')
        Esort=np.multiply(self.eigenval,hartreetocm)
        if saveeigen and not self.petsc:
            eigfile='tmp.eig.h5'
            from os.path import isfile
            if isfile(eigfile):
                if not self.saveeigrequested:
                    if 'y' not in input('outfile (tmp.eig) exists, overwrite? [y,N]').lower():
                        saveeigen=False
                    self.saveeigrequested=True
            if saveeigen:
                if np.abs(np.subtract(mw1,mw2))>1.0E-07 or self.mingrid:
                    gridx, gridy, pot = np.ravel(grid_x),np.ravel(grid_y),np.ravel(vfit)
                else:
                    gridx, gridy = np.array(r[0],dtype=eval(numpy_precision)), np.array(r[1],dtype=eval(numpy_precision)) 
                    pot= np.array(Energies,dtype=eval(numpy_precision)) 
                import h5py
                f = h5py.File(eigfile,'w')
                f.create_dataset('x',data=gridx)
                f.create_dataset('y',data=gridy)
                f.create_dataset('z',data=pot)
                f.create_dataset('eigenvec',data=eigenvec)
                f.create_dataset('eigenval',data=self.eigenval)
                f.close()
        if self.plotit and not self.petsc:
            if np.abs(np.subtract(mw1,mw2))>1.0E-07 or self.mingrid:
                plot2dgrid(grid_x,grid_y,vfit,wavenumber=True,title='Potential Energy Contours')
                eigenvectoplot=(int(input('number of eigenvectors to plot:')))
                if eigenvectoplot>0:
                    for i in range(eigenvectoplot):
                        if i==eigenvectoplot-1:
                            plot2d(np.ndarray.flatten(grid_x),np.ndarray.flatten(grid_y),np.square(eigenvec[i]),\
                                    title='eigenvec {0} with energy {1:.3f}'.format(i,Esort[i]),block=True)
                        else:
                            plot2d(np.ndarray.flatten(grid_x),np.ndarray.flatten(grid_y),np.square(eigenvec[i]),\
                                    title='eigenvec {0} with energy {1:.3f}'.format(i,Esort[i]),block=False)
            else:
                plot2d(r[0],r[1],Energies,wavenumber=True,title='Potential Energy Contours')
                eigenvectoplot=(int(input('number of eigenvectors to plot:')))
                if eigenvectoplot>0:
                    for i in range(eigenvectoplot):
                        if i==eigenvectoplot-1:
                            plot2d(r[0],r[1],np.square(eigenvec[i]),title='eigenvec {0} with energy {1:.3f}'.format(i,Esort[i]),block=True)
                        else:
                            plot2d(r[0],r[1],np.square(eigenvec[i]),title='eigenvec {0} with energy {1:.3f}'.format(i,Esort[i]))

def loadeigen(eigfile='tmp.eig.h5',eigenvectoplot=1):
    import numpy as np
    import h5py
    eigbase=eigfile.split('.')[0]
    data=h5py.File(eigfile,'r')
    x=data['x'][:]
    y=data['y'][:]
    z=data['z'][:]
    eigenvec=data['eigenvec'][:]
    eigenval=data['eigenval'][:]
    Esort=np.multiply(eigenval,hartreetocm)
    plot2d(x,y,z,wavenumber=True,title='Potential Energy Contours',save=eigbase+'.pot.pdf')
    if eigenvectoplot>0:
        for i in range(eigenvectoplot):
            if i==eigenvectoplot-1:
                plot2d(x,y,eigenvec[i],title='eigenvec {0} with energy {1:.3f}'.format(i,Esort[i]),save=eigbase+'.eig.'+str(i)+'.pdf',includegrid=False)
            else:
                plot2d(x,y,eigenvec[i],title='eigenvec {0} with energy {1:.3f}'.format(i,Esort[i]),save=eigbase+'.eig.'+str(i)+'.pdf',includegrid=False)

def plot2dgrid(x,y,z,wavenumber=False,angular=False,norm=False,block=False,legend=True,title='2d filled contour plot',save='tmp.pdf'):
    """ 2d grids in with their corresponding z, e.g. np.shape (x_dim_len,y_dim_len) for x, y, and z"""
    import numpy as np
    import matplotlib.pyplot as plt
    if norm:
        z=np.divide(z,np.subtract(np.amax(z),np.amin(z)))
    if wavenumber:
        z=np.subtract(z,np.min(z))
        z=z*219474.6313717 
        levs=range(0,10000,100)
    else:
        levs=np.linspace(float(np.amin(z)),float(np.amax(z)),100)
    plt.figure()
    cp=plt.contourf(x,y,z,cmap=(plt.cm.gnuplot),origin='lower',levels=levs)
    plt.title(title)
    if legend:
        CB = plt.colorbar(cp, shrink=0.8, extend='both')
        l, b, w, h = plt.gca().get_position().bounds
        ll, bb, ww, hh = CB.ax.get_position().bounds
        CB.ax.set_position([ll, b + 0.1*h, ww, h*0.8])
    if save=='tmp.pdf':
        plt.show(block=block)
    else:
        plt.savefig(save)

def plot2d(x,y,z,wavenumber=False,wavenumbercutoff=10000,angular=True,norm=False,block=False,legend=True,includegrid=False,title='2d filled contour plot',save='tmp.pdf'):
    """ Flat x,y,z as input"""
    import numpy as np
    import matplotlib.mlab as ml
    import matplotlib.pyplot as plt
    if norm and wavenumber:
        raise ValueError("normalization and wavenumber plotting muturally inconsistent ")
    if norm:
        z=np.divide(z,np.subtract(np.max(z),np.min(z)))
    if wavenumber:
        z=np.subtract(z,np.min(z))
        z=z*219474.6313717 
        levs=range(0,wavenumbercutoff,int(min(100,wavenumbercutoff/50)))
#        levs=range(0,wavenumbercutoff,100)
    else:
        levs=np.linspace(float(np.min(z)),float(np.max(z)),100)
    xlen=len(set(x))
    ylen=len(set(y))
    if angular:
        x=(x*180/np.pi)
        y=(y*180/np.pi)
    xi=np.linspace(min(x),max(x),xlen)
    yi=np.linspace(min(y),max(y),ylen)
    zi=ml.griddata(x,y,z,xi,yi,interp='nn')
    plt.figure()
#    cp=plt.contour(xi,yi,zi,cmap=(plt.cm.gnuplot),origin='lower',levels=levs)
    cp=plt.contourf(xi,yi,zi,cmap=(plt.cm.gnuplot),origin='lower',levels=levs)
    if includegrid:
        plt.plot(x,y,linestyle='none',marker='o')
    plt.title(title)
    if legend:
        CB = plt.colorbar(cp, shrink=0.8, extend='both')
#        CB = plt.colorbar(cp, orientation='horizontal', shrink=0.8)
        l, b, w, h = plt.gca().get_position().bounds
        ll, bb, ww, hh = CB.ax.get_position().bounds
        CB.ax.set_position([ll, b + 0.1*h, ww, h*0.8])
    if save=='tmp.pdf':
        plt.show(block=block)
    else:
        plt.savefig(save)
        plt.close()

def H_array_petsc(pts=5,coordtype=['r'],mass=[0.5],qmin=[1.0],qmax=[2.0],V=[],numeig=6,printpetsc=False):
    import sys, slepc4py
    slepc4py.init(sys.argv)
    from petsc4py import PETSc
    from slepc4py import SLEPc
    totpts=1
    for x in range(len(pts)):
        totpts=totpts*pts[x]
    opts=PETSc.Options()
    from mpi4py import MPI
    import numpy as np
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    totalproc= comm.size
    A=PETSc.Mat().create()
    A.setType('sbaij')
#    A.setType('mpisbaij')
    A.setSizes([(None,totpts),(None,totpts)])
    A.setFromOptions()
    A.setUp()
    diagvec=PETSc.Vec().createSeq(V.shape[0])
    diagvec.setValues(range(V.shape[0]),V)
    diagvec.assemble()
    rstart,rend=A.getOwnershipRange()
    ncoord=len(coordtype)
    if ncoord==1:
        qmin=[qmin]
        qmax=[qmax]
        D1, mw1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0],coordtype=coordtype[0])
        if np.less(pts[0],255):
            inttype='np.uint8'
        else:
            inttype='np.uint16'
        indices=np.array(np.meshgrid(np.arange(pts[0]),np.arange(pts[0])),dtype=eval(inttype)).T.reshape(-1,2)
    if ncoord==2:
        D1, mw1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0],coordtype=coordtype[0])
        D2, mw2=H_array_1d(pts[1],mass=mass[1],qmin=qmin[1],qmax=qmax[1],coordtype=coordtype[1])
        if np.abs(np.subtract(mw1,mw2))>1.0E-07:
            from sys import exit
            print('mass weighted coordinate spacing unequal as specified, stopping DVR')
            exit('mass weighted coordinate spacing unequal as specified, stopping DVR')
        if np.less(pts[0],255) and np.less(pts[1],255):
            inttype='np.uint8'
        else:
            inttype='np.uint16'
        indices=np.array(np.meshgrid(np.arange(pts[0]),np.arange(pts[1]),np.arange(pts[1]),np.arange(pts[0])),dtype=eval(inttype)).T.reshape(-1,4)
    if ncoord==3:
        raise ValueError("not implemented for %d dimensions" % (ncoord))
        D1, mw1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0],coordtype=coordtype[0])
        D2, mw2=H_array_1d(pts[1],mass=mass[1],qmin=qmin[1],qmax=qmax[1],coordtype=coordtype[1])
        D3, mw3=H_array_1d(pts[2],mass=mass[2],qmin=qmin[2],qmax=qmax[2],coordtype=coordtype[2])
    indices=np.array(np.split(indices,ncoord*2,axis=1),dtype=eval(inttype))
    k=0
    if ncoord==1:
        for i in range(rstart,rend):
            for j in range(rstart,rend):
                if i==j:
                    A[i,i]=np.add(D1[i,j],A[i,i])
                else:
                    A[i,j]=D1[i,j]
    elif ncoord==2:
        D1add=np.equal(indices[0,:],indices[3,:])
        D2add=np.equal(indices[1,:],indices[2,:])
        ijindex=np.array(np.meshgrid(np.arange(totpts),np.arange(totpts)),dtype=np.uint64).T.reshape(-1,2)
        ijrange=np.squeeze(np.where(np.logical_or(D1add,D2add))[0])
        localijrange=np.array_split(ijrange,totalproc)[rank]
        if rank==0 and printpetsc:
            print('Sparse % {0}'.format(np.multiply(100,np.divide(np.add(V.shape[0],ijrange.shape[0]),np.square(totpts)))))
        for x in np.nditer(localijrange):
            if D1add[x]: 
                if D2add[x]:
                        A[ijindex[x,0],ijindex[x,1]]=np.add(np.add(D1[indices[0,x],indices[3,x]], D2[indices[1,x],indices[2,x]]),diagvec[ijindex[x,0]])
                else:
                    A[ijindex[x,0],ijindex[x,1]]=D2[indices[1,x],indices[2,x]]
            else:
                A[ijindex[x,0],ijindex[x,1]]=D1[indices[0,x],indices[3,x]]
    else:
        raise ValueError("not implemented for %d dimensions" % (ncoord))
    A.assemble()
    E = SLEPc.EPS(); E.create()
    E.setOperators(A)
    E.setProblemType(SLEPc.EPS.ProblemType.HEP)
    E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)
    E.setDimensions(numeig)
    E.setFromOptions()
    E.solve()
    if rank==0 and printpetsc:
        E.view()
        its = E.getIterationNumber()
        Print = PETSc.Sys.Print
        Print("Number of iterations of the method: %d" % its)
    if rank==0 and printpetsc:
        eps_type = E.getType()
        Print("Solution method: %s" % eps_type)
    nev, ncv, mpd = E.getDimensions()
    tol, maxit = E.getTolerances()
    nconv = E.getConverged()
    eigenval=np.zeros(nconv,dtype=np.float64)
    if nconv > 0:
        vr, wr = A.getVecs()
        vi, wi = A.getVecs()
        for i in range(nconv):
            k = E.getEigenpair(i, vr, vi)
            error = E.computeError(i)
            eigenval[i]=k.real
    return (eigenval)

def H_array(pts=5,coordtype=['r'],mass=[0.5],qmin=[1.0],qmax=[2.0],V=[]):
    """ input 
    pts=points , coordtype (dict values r, phi, theta), qmin=list(qminimia), qmax=list(qmaxima), V=potential
    output
    2d DVR array based on
    Kinetic Energy Array (dimensionality=2): see Eq A6a and A6b of JCP 96, 1982 (1992): note 
    constants are defined in constants module globally earlier
    The Hamiltonian has been converted to atomic units, e.g.
    H =  - [1/(2 am)] d^2/dx^2 + v(x)
    ncoord must be passed, simplest is 1 for a 1D potential
    note a and b always occur outside of the potential!!!
    for 0 to 2 pi only one of the two is included
    The massweighted spacing. All coordinates must use the same
    pts is the number of points per coordinate
    mass given in amu; converted to atomic units here
    For angular coordinates this is instead  a moment of intertia in units of amu * bohr^2 """
    import numpy as np
    np.set_printoptions(suppress=False,threshold=np.nan,linewidth=np.nan)
    ncoord=len(coordtype)
    totpts=1
    for x in range(len(pts)):
        totpts=totpts*pts[x]
    A=np.zeros((totpts,totpts),dtype=eval(numpy_precision))
    if ncoord==1:
        qmin=[qmin]
        qmax=[qmax]
        D1, mw1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0],coordtype=coordtype[0])
        if np.less(pts[0],255):
            inttype='np.uint8'
        else:
            inttype='np.uint16'
        indices=np.array(np.meshgrid(np.arange(pts[0]),np.arange(pts[0])),dtype=eval(inttype)).T.reshape(-1,2)
    if ncoord==2:
        D1, mw1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0],coordtype=coordtype[0])
        D2, mw2=H_array_1d(pts[1],mass=mass[1],qmin=qmin[1],qmax=qmax[1],coordtype=coordtype[1])
        if np.abs(np.subtract(mw1,mw2))>1.0E-07:
            from sys import exit
            print('mass weighted coordinate spacing unequal as specified, stopping DVR')
            exit('mass weighted coordinate spacing unequal as specified, stopping DVR')
        if np.less(pts[0],255) and np.less(pts[1],255):
            inttype='np.uint8'
        else:
            inttype='np.uint16'
        indices=np.array(np.meshgrid(np.arange(pts[0]),np.arange(pts[1]),np.arange(pts[1]),np.arange(pts[0])),dtype=eval(inttype)).T.reshape(-1,4)
    if ncoord==3:
        raise ValueError("not implemented for %d dimensions" % (ncoord))
        D1, mw1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0],coordtype=coordtype[0])
        D2, mw2=H_array_1d(pts[1],mass=mass[1],qmin=qmin[1],qmax=qmax[1],coordtype=coordtype[1])
        D3, mw3=H_array_1d(pts[2],mass=mass[2],qmin=qmin[2],qmax=qmax[2],coordtype=coordtype[2])
#        indices=np.array(np.meshgrid(np.arange(ptspercoord),np.arange(ptspercoord),np.arange(ptspercoord),np.arange(ptspercoord)\
#                ,np.arange(ptspercoord),np.arange(ptspercoord))).T.reshape(-1,6)
    indices=np.array(np.split(indices,ncoord*2,axis=1),dtype=eval(inttype))
    it= np.nditer(A, flags=['c_index'], op_flags=['writeonly'])
    k=0
    if ncoord==1:
        while not it.finished:
            i,i1 =indices[0,it.index], indices[1,it.index]
            if i==i1:
                it[0]=np.add(D1[i,i1],V[i])
            else:
                it[0]=D1[i,i1]
            it.iternext()
    elif ncoord==2:
        D1add=np.equal(indices[0,:],indices[3,:])
        D2add=np.equal(indices[1,:],indices[2,:])
        ijindex=np.array(np.meshgrid(np.arange(totpts),np.arange(totpts)),dtype=np.uint64).T.reshape(-1,2)
        ijrange=np.squeeze(np.where(np.logical_or(D1add,D2add))[0])
        for x in np.nditer(ijrange):
            if D1add[x]: 
                if D2add[x]:
                        A[ijindex[x,0],ijindex[x,1]]=np.add(np.add(D1[indices[0,x],indices[3,x]], D2[indices[1,x],indices[2,x]]),V[ijindex[x,0]])
                else:
                    A[ijindex[x,0],ijindex[x,1]]=D2[indices[1,x],indices[2,x]]
            else:
                A[ijindex[x,0],ijindex[x,1]]=D1[indices[0,x],indices[3,x]]
    else:
        raise ValueError("not implemented for %d dimensions" % (ncoord))
    return A

def H_array_1d(pts=5,coordtype='r',mass=0.5,qmin=1.0,qmax=2.0):
    import numpy as np
    n=pts
    A=np.zeros((n,n),dtype=eval(numpy_precision))
# In atomic units
# One has been added to i and j inside to make this consistent with paper
    mass_conv=mass*(amu/e_mass)
    rlen=0
    for i in range(pts):
        for j in range(pts):
            if coordtype=='r':
                n1=pts+1
                """ rlen here is (max-min) of potential, scaled to b-a by adding two more points!"""
                rlen=np.multiply(np.subtract(qmax,qmin),np.divide(np.add(float(pts),1.0),np.subtract(float(pts),1.0)))
                prefactor=np.divide(np.power(np.pi,2),np.multiply(np.multiply(4,mass_conv),np.power(rlen,2)))
                if i==j:
                    A[i,i]=np.multiply(prefactor,np.subtract(\
                            np.divide(np.add(np.multiply(2,np.power(n1,2)),1),3),np.power(np.sin(np.divide(np.multiply(np.add(i,1),np.pi),n1)),-2)))
                else:
                    A[i,j]=np.multiply(np.multiply(prefactor,np.power(-1,np.subtract(i,j)))\
                            ,np.subtract(np.power(np.sin(np.divide(np.multiply(np.pi,np.subtract(i,j)) , np.multiply(2 ,n1) )),-2) , \
                            np.power(np.sin(np.divide(np.multiply(np.pi,np.add(i,np.add(j,2))) , np.multiply(2 , n1))),-2)))
# 0 to 2pi in appendix A section 4
            elif coordtype=='phi':
                prefactor=(1.0)/(2*mass_conv)
                m=int(np.divide(pts,2))
                if (2*m+1)!=pts:
                    from sys import exit
                    exit('in phi coordinate 2m+1 != n, must use odd number of points')
                if i==j:
                    A[np.sum(i),np.sum(j)]=np.multiply(prefactor,np.divide(np.multiply(m,np.add(m,1)),3))
                else:
                    cosij=np.cos(np.divide(np.multiply(np.pi,np.subtract(i,j)),n))
                    A[np.sum(i),np.sum(j)]=\
                            np.multiply(np.multiply(np.power(-1,np.subtract(i,j)),prefactor),np.divide(cosij,np.multiply(2,np.subtract(1,np.power(cosij,2)))))
            elif coordtype=='theta':
                n1=pts+1
                prefactor=np.divide(1.0,np.multiply(4,mass_conv))
                if i==j:
                    A[i,i]=np.multiply(prefactor,np.subtract(\
                            np.divide(np.add(np.multiply(2,np.power(n1,2)),1),3),np.power(np.sin(np.divide(np.multiply(np.add(i,1),np.pi),n1)),-2)))
                else:
                    A[i,j]=np.multiply(np.multiply(prefactor,np.power(-1,np.subtract(i,j)))\
                            ,np.subtract(np.power(np.sin(np.divide(np.multiply(np.pi,np.subtract(i,j)) , np.multiply(2 ,n1) )),-2) , \
                            np.power(np.sin(np.divide(np.multiply(np.pi,np.add(i,np.add(j,2))) , np.multiply(2 , n1))),-2)))
                pass
            else:
                from sys import exit
                exit('coordinate type not recongized')
    mws=mwspace(rlen=rlen,mass=mass,coordtype=coordtype,pts=pts)
#    if coordtype=='r' or coordtype=='theta':
#        print('{1} has prefactor of {0:.4e} with b-a={2}'.format(prefactor, coordtype,rlen))
#    else:
#        print('{1} has prefactor of {0:.4e}'.format(prefactor/8, coordtype))
    return (A,mws)

def openandread(filename):
    """ Opens a file and returns the lines. Upon UnicodeDecodeError will try latin1 encoding and return those lines. Encodes all lines as utf8 and strips linebreaks before returning"""
    import sys,os.path
    if os.path.isfile(filename):
        txt=open(filename,'r')
    else:
        sys.exit('file not found')
    try:
        lines=txt.readlines()
    except UnicodeDecodeError:
        import codecs
        txt = codecs.open(filename,'r',encoding='latin1')
        lines=txt.readlines()
    for x in range(len(lines)):
        lines[x]=lines[x].encode('utf8').decode('utf8').strip()
    return lines

def jacobi(A,b,N=25,x=None):
    """Solves the equation Ax=b via the Jacobi iterative method."""
    from numpy import array, zeros, diag, diagflat, dot
    # Create an initial guess if needed
    if x is None:
        x = zeros(len(A[0]))
    # Create a vector of the diagonal elements of A
    # and subtract them from A
    D = diag(A)
    R = A - diagflat(D)
    # Iterate for N times
    for i in range(N):
        x = (b - dot(R,x)) / D
    return x

def return1dsplinemin(x,y):
    """Returns the minimum from a spline function"""
    from scipy.interpolate import splrep, splder, sproot, splev
    spline=splrep(x, y, s=0,k=4)
    xmin=sproot(splder(spline) )
    return (xmin, splev(xmin,spline,der=0)) 
#    return spline.interpolate.derivative().roots()

def return2dspline(c,spline,dx,dy):
    """ Returns value of 2d spline, for requested derivitive"""
    return spline.ev(c[0],c[1],dx=dx,dy=dy)

def return2dsplinetotder(c,spline,dx,dy):
    """ Returns value of 2d spline, for sum of dx=1 and dy=1. Ignores dx dy args (because of scipy.minimize behavior for jacobian"""
    from numpy import array
    return array(([spline.ev(c[0],c[1],dx=1,dy=0),spline.ev(c[0],c[1],dx=0,dy=1)]))

def func1d(x,c,a1,b1,a2,b2,a3,b3,a4,b4):
    from numpy import cos, sin, add, multiply, power
    return add(c, add( add( add( \
            add(multiply(cos(x),a1),multiply(sin(x),b1)) , \
            add(multiply(cos(multiply(x,2)),a2),multiply(sin(multiply(x,2)),b2))),\
            add(multiply(cos(multiply(x,3)),a3),multiply(sin(multiply(x,3)),b3))),\
            add(multiply(cos(multiply(x,4)),a4),multiply(sin(multiply(x,4)),b4))),\
            )

def func1dder1(x,c,a1,b1,a2,b2,a3,b3,a4,b4):
    from numpy import cos, sin, add, multiply
    return  add( add( add(\
            add(multiply(-sin(x),a1),multiply(cos(x),b1)) , \
            add(multiply(-sin(multiply(x,2)),multiply(a2,2)),multiply(cos(multiply(x,2)),multiply(b2,2)))),\
            add(multiply(-sin(multiply(x,3)),multiply(a3,3)),multiply(cos(multiply(x,3)),multiply(b3,3)))),\
            add(multiply(-sin(multiply(x,4)),multiply(a4,4)),multiply(cos(multiply(x,4)),multiply(b4,4))),\
            )

def func1dder2(x,c,a1,b1,a2,b2,a3,b3,a4,b4):
    from numpy import cos, sin, add, multiply
    return  add( add( add(\
            add(multiply(-cos(x),a1),multiply(cos(x),b1)) , \
            add(multiply(-cos(multiply(x,2)),multiply(a2,4)),multiply(-sin(multiply(x,2)),multiply(b2,4)))),\
            add(multiply(-cos(multiply(x,3)),multiply(a3,9)),multiply(-sin(multiply(x,3)),multiply(b3,9)))),\
            add(multiply(-cos(multiply(x,4)),multiply(a4,16)),multiply(-sin(multiply(x,4)),multiply(b4,16))),\
            )

def plot1dfit(x0,y0,xgrid,yfit,block=True,legend=['Points','Fit'],title='Fit Potential'):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(x0,y0,linestyle='none',marker='o')
    plt.plot(xgrid,yfit,linestyle='solid',marker='None')
    plt.axis()
    plt.legend(legend)
    plt.title(title)
    plt.show(block=block)

    
def mwspace(coordtype='r',rlen=1.0,mass=1.0,pts=2):
    from numpy import multiply, power, divide, pi , add
    if coordtype=='r':
        mwspace=multiply(power(divide(rlen,add(pts,1)),2),mass)
    elif coordtype=='theta':
        mwspace=multiply(power(divide(pi,add(pts,1)),2),mass)
    elif coordtype=='phi':
        mwspace=multiply(power(divide(multiply(2,pi),pts),2),mass)
    return mwspace

def main():
    constants(CODATA_year=2010)
    import numpy as np
    import sys

    """ overloaded call... I may split this out later"""
    if len(sys.argv)>1 and 'h5' in sys.argv[1]:
        if len(sys.argv)>2:
            loadeigen(eigfile=sys.argv[1],eigenvectoplot=int(sys.argv[2]))
        else:
            loadeigen(eigfile=sys.argv[1],eigenvectoplot=int(input('number of eigenvectors to plot:'))) 
    else:
        pot=potential()
        if len(sys.argv)>1:
            pot.readpotential(inp=sys.argv[1])
        else:
            pot.readpotential(inp=input('Give the file with the potential: '))
#        pot.xlsx()
#       pot.fit1dpot()
        pot.solve()
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if len(pot.fiteignum)>0 and rank==0:
            pot.fitfundamental()
        if pot.printeigenval==True and rank==0:
            pot.printeigenvals()

# jacobian stuff
#    A = array([[2.0,1.0],[5.0,7.0]])
#    b = array([11.0,13.0])
#    guess = array([1.0,1.0])
#    import numpy as np
#    sol = jacobi(A,b,N=25,x=guess)
#    from pprint import pprint
#    pprint(sol)
#    from sys import exit
#    exit()

if __name__=="__main__":
    main()
