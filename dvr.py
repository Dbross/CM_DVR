#!/usr/bin/env python3
""" This program was written following 
Citation: The Journal of Chemical Physics 96, 1982 (1992); doi: 10.1063/1.462100
as a guide
Can use  scipy.sparse.linalg.eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')[source] as a way of doing Lanczos...."""
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (bytes, str, open, super, range, zip, round, input, int, pow, object)

def constants(CODATA_year=2010):
    """ CODATA constants used in program and other global defintions including number of points and numpy datatype for reals"""
    global numpy_precision, num_points, num_print_digits, plotit
    numpy_precision="np.float64"
    num_points=101
    num_print_digits=3
    plotit=True
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

def readpotential(inp,r_units='bohr'):
#potential should have coordinate and units as main input
# anyline starting with ! or # is commented out
    commentoutchars=['!','#']
    lines=openandread(inp)
    r=[]
    energy=[]
    coordtypes=[]
    types=['angular_2pi','radial','angular_pi']
    mass=0
    if r_units.lower()=='bohr':
        r_unitconversion=1.0
    elif r_units.lower()=='angstrom':
        r_unitconversion=(1.0/bohr)
#        r_unitconversion=1.88972613
    else:
        from sys import exit
        exit('No valid units given for length')
    import re
    numeric_const_pattern = r"""[-+]?(?: (?: \d* \. \d+ ) | (?: \d+ \.? ) ) (?: [Ee] [+-]? \d+ ) ?"""
    rx=re.compile(numeric_const_pattern,re.VERBOSE)
    emin=0
    for x in lines:
        if x[0] not in commentoutchars:
            if 'mass' in x.lower() and mass==0:
                mass=rx.findall(x)
                for x in range(len(mass)):
                    mass[x]=float(mass[x])
                print('using reduced mass of {0} amu.'.format(mass))
            elif 'mass' in x.lower():
                from sys import exit
                exit('mass defined in potential twice')
            elif 'bohr' in x.lower():
                print('reading potential as bohr, this should only be set once.')
                r_unitconversion=1.0
            elif 'emin' in x.lower():
                emin=rx.findall(x)
                if len(emin)>1:
                    for x in emin:
                        x=float(x)
                    print('found multiple energy minimum {0} using {1} as minimum'.format(emin,min(emin)))
                    emin=min(emin)
                else:
                    emin=float(emin[0])
                print('shifting energy minimum to {} a.u.'.format(emin))
            elif types[0] in x.lower() or types[1] in x.lower() or types[2] in x.lower():
                typelist=x.lower().replace(',',' ').split()
                for y in typelist:
                    if y in types:
                        coordtypes.append(y)
                #if len(coordtypes)<len(coords):
                #    from sys import exit
                #    exit('More coordinate types than coordinates')
            elif 'angstrom' in x.lower():
                print('reading potential as angstrom, this should only be set once.')
                r_unitconversion=(1.0/bohr)
            else:
                linesplit=x.split()
                rtmp=[]
                for x in range(0,len(linesplit)-1):
                    rtmp.append(float(linesplit[x])*r_unitconversion)
                r.append(rtmp)
                energy.append(float(linesplit[len(linesplit)-1]))
                if len(r[-1])!=len(coordtypes):
                    from sys import exit
                    print(x)
                    exit('number of coordinates given inconsistent with coordinate types given')
        else:
            print('{0} commented out'.format(x[1:].strip()))
    rnew=[]
    for i in range(len(r[0])):
        rtmp=[r[j][i] for j in range(len(r))]
        rnew.append(rtmp)
    if len(mass)!=len(coordtypes):
        print('{0} masses given and {1} coordinate types give'.format(len(mass),len(coordtypes)))
        from sys import exit
        exit()
    if len(mass)!=len(coordtypes):
        print('{0} masses given and {1} coordinate types give'.format(len(mass),len(coordtypes)))
        from sys import exit
        exit()
    for x in range(len(energy)):
        energy[x]=float(energy[x])-emin
    return (rnew,energy, mass, coordtypes)


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

def cubicspline(x,y):
    """ Performs a cubic spline with no smoothing and returns the spline function"""
    from scipy.interpolate import splrep
    return splrep(x, y, s=0)

def returnsplinevalue(spline,xnew):
    """ returns the value(s) of a spline function at a given x(s)"""
    from scipy.interpolate import splev
    return splev(xnew, spline, der=0)

def returnsplinemin(x,y):
    """Returns the minimum from a spline function"""
    from scipy.interpolate import splrep, splder, sproot
    spline=splrep(x, y, s=0,k=4)
    xmin=sproot(splder(spline) )
    return (xmin,returnsplinevalue(spline,xmin))
#    return spline.interpolate.derivative().roots()

def H_array_1d(pts=5,coordtype='r',mass=0.5,qmin=1.0,qmax=2.0):
    import numpy as np
    n=pts
    A=np.zeros((n,n),dtype=eval(numpy_precision))
# In atomic units
# One has been added to i and j inside to make this consistent with paper
    mass_conv=mass*(amu/e_mass)
    for i in range(pts):
        for j in range(pts):
            if coordtype=='r':
                n1=pts+1
                dq=np.multiply(np.subtract(qmax,qmin),np.divide(np.add(float(pts),1.0),np.subtract(float(pts),1.0)))
                prefactor=np.divide(np.power(np.pi,2),np.multiply(np.multiply(4,mass_conv),np.power(dq,2)))
                if i==j:
                    A[i,i]=np.multiply(prefactor,np.subtract(\
                            np.divide(np.add(np.multiply(2,np.power(n1,2)),1),3),np.power(np.sin(np.divide(np.multiply(np.add(i,1),np.pi),n1)),-2)))
                else:
                    A[i,j]=np.multiply(np.multiply(prefactor,np.power(-1,np.subtract(i,j)))\
                            ,np.subtract(np.power(np.sin(np.divide(np.multiply(np.pi,np.subtract(i,j)) , np.multiply(2 ,n1) )),-2) , \
                            np.power(np.sin(np.divide(np.multiply(np.pi,np.add(i,np.add(j,2))) , np.multiply(2 , n1))),-2)))
# 0 to 2pi in appendix A section 4
            elif coordtype[x]=='phi':
                prefactor=(1.0)/(2*mass_conv)
                m=int(np.divide(pts,2))
                if (2*m+1)!=pts:
                    from sys import exit
                    exit('in phi coordinate 2m+1 != n, must use odd number of points')
                if i==j:
                    A[np.sum(i,gridstart),np.sum(j,gridstart)]=np.add(np.multiply(prefactor,np.divide(np.multiply(m,np.add(m,1)),3)),V[i])
                else:
                    cosij=np.cos(np.divide(np.multiply(np.pi,np.subtract(i,j)),n))
                    A[np.sum(i,gridstart),np.sum(j,gridstart)]=\
                            np.multiply(np.multiply(np.power(-1,np.subtract(i,j)),prefactor),np.divide(cosij,np.multiply(2,np.subtract(1,np.power(cosij,2)))))
            elif coordtype[x]=='theta':
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
    return A

def H_array(pts=5,coordtype=['r'],mass=[0.5],dq=0.001,qmin=[1.0],qmax=[2.0],V=[]):
    """ Kinetic Energy Array (dimensionality=2): see Eq A6a and A6b of JCP 96, 1982 (1992): note 
    constants are defined in constants module globally earlier
    The Hamiltonian has been converted to atomic units, e.g.
    H =  - [1/(2 am)] d^2/dx^2 + v(x)
    ncoord must be passed, simplest is 1 for a 1D potential
    dq is the massweighted spacing. All coordinates must use the same
    pts is the number of points per coordinate
    mass given in amu; converted to atomic units here"""
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
        D1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0])
        indices=np.array(np.meshgrid(np.arange(pts[0]),np.arange(pts[0]))).T.reshape(-1,2)
    if ncoord==2:
        D1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0])
        D2=H_array_1d(pts[1],mass=mass[1],qmin=qmin[1],qmax=qmax[1])
        indices=np.array(np.meshgrid(np.arange(pts[0]),np.arange(pts[1]),np.arange(pts[1]),np.arange(pts[0]))).T.reshape(-1,4)
    if ncoord==3:
        D1=H_array_1d(pts[0],mass=mass[0],qmin=qmin[0],qmax=qmax[0])
        D2=H_array_1d(pts[1],mass=mass[1],qmin=qmin[1],qmax=qmax[1])
        D3=H_array_1d(pts[2],mass=mass[2],qmin=qmin[2],qmax=qmax[2])
#        indices=np.array(np.meshgrid(np.arange(ptspercoord),np.arange(ptspercoord),np.arange(ptspercoord),np.arange(ptspercoord)\
#                ,np.arange(ptspercoord),np.arange(ptspercoord))).T.reshape(-1,6)
#    from itertools import product
#    indices=np.array([item for item in product(range(ptspercoord),repeat=ncoord*2)],dtype=int)
    indices=np.split(indices,ncoord*2,axis=1) 
    it= np.nditer(A, flags=['c_index'], op_flags=['writeonly'])
    k=0
    if ncoord==1:
        while not it.finished:
            i=indices[0][it.index]
            i1=indices[1][it.index]
            if i==i1:
                it[0]=np.add(D1[i,i1],V[i])
            else:
                it[0]=D1[i,i1]
            it.iternext()
    elif ncoord==2:
        while not it.finished:
            x=it.index
            i=indices[0][x]
            i1=indices[3][x]
            j=indices[1][x]
            j1=indices[2][x]
            if i==i1 and j==j1:
                it[0]=np.add(np.add(D1[i,i1],D2[j,j1]),V[k])
                k+=1
            elif i==i1:
                it[0]=D2[j,j1]
            elif j==j1:
                it[0]=D1[i,i1]
            it.iternext()
    return A

def spline1dpot(pts,mass,coordtypes,Energies_raw,r_raw):
    import numpy as np
    xmin,emin=returnsplinemin(r_raw[0],Energies_raw)
    r=r_raw-np.min(xmin)
    Energies=Energies_raw-np.min(emin)
    Ener_spline=cubicspline(r[0],Energies)
    xnew = np.linspace(min(r[0]),max(r[0]), num=num_points)
    vfit=returnsplinevalue(Ener_spline,xnew)
    pts=[]
    for x in range(len(r)):
        pts.append(len(np.unique(r[x])))
    Ham=H_array(pts=pts,mass=mass,V=Energies,qmax=np.amax(r,axis=1),qmin=np.amin(r,axis=1),coordtype=coordtypes)
    eigenval, eigenvec=np.linalg.eig(Ham)
    Esort=np.sort(eigenval*hartreetocm)
# plotting stuff 
    Etoprint=int(len(Esort)/2)
    if plotit:
        vfitcm=vfit*hartreetocm
        import matplotlib.pyplot as plt
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
#    plt.legend(['Points', 'Cubic Spline'])
        plt.title('Cubic-spline interpolation')
        plt.axis()
        plt.show()
    for x in range(Etoprint):
        print('{0:.{1}f}'.format(round(Esort[x],num_print_digits),num_print_digits))
    return eigenval 

def main():
    constants(CODATA_year=2010)
    import numpy as np
    import sys
    if len(sys.argv)>1:
        potential=readpotential(sys.argv[1],r_units='angstrom')
    else:
        potential=readpotential(input('Give the file with the potential: '),r_units='angstrom')
    r=np.array(potential[0],dtype=eval(numpy_precision))
    Energies=np.array(potential[1],dtype=eval(numpy_precision))
    mass=potential[2]
    coordtypes=potential[3]
    """ Radial terminates with PiB walls; angular_2pi repeats; angular_pi terminates with PiB walls at 0 and pi"""
    coordtypedict={'radial': 'r', 'angular_2pi': 'phi', 'angular_pi': 'theta'}
    for x in range(len(coordtypes)):
        coordtypes[x]=coordtypedict[coordtypes[x]]
    pts=[]
    for x in range(len(r)):
        pts.append(len(np.unique(r[x])))
    if len(coordtypes)==1:
        eigenval= spline1dpot(pts,mass,coordtypes,Energies,r)
#    elif len(coordtypes)==2:
#        eigenval= spline2dpot(pts,mass,coordtypes,Energies,r)
    else:
        Ham=H_array(pts=pts,mass=mass,V=Energies,qmax=np.amax(r,axis=1),qmin=np.amin(r,axis=1),coordtype=coordtypes)
        eigenval, eigenvec=np.linalg.eig(Ham)
        Esort=np.sort(eigenval*hartreetocm)
        Etoprint=int(len(Esort)/2)
        for x in range(Etoprint):
            print('{0:.{1}f}'.format(round(Esort[x],num_print_digits),num_print_digits))

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
