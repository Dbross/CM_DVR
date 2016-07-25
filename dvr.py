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
    for x in lines:
        if x[0] not in commentoutchars:
            if 'mass' in x.lower() and mass==0:
                mass=re.findall("\d+\.\d+",x)
                if len(mass)>1:
                    from sys import exit
                    exit('mass defined in potential twice')
                print('using reduced mass of {0} amu.'.format(mass[0]))
                mass=float(mass[0])
            elif 'mass' in x.lower():
                from sys import exit
                exit('mass defined in potential twice')
            elif 'bohr' in x.lower():
                print('reading potential as bohr, this should only be set once and the line it is on is ignored.')
                r_unitconversion=1.0
            elif 'angstrom' in x.lower():
                print('reading potential as angstrom, this should only be set once and the line it is on is ignored.')
                r_unitconversion=(1.0/bohr)
            else:
                linesplit=x.split()
                if len(linesplit)==2:
                    r.append(float(linesplit[0])*r_unitconversion)
                    energy.append(float(linesplit[1]))
                else:
                    from sys import exit
                    print(x)
                    exit('line in potential not properly parsed')
    return (r,energy, mass)


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

def H_array(ncoord=1,pts=5,coordtype=['r'],mass=0.5,dq=0.001,qmax=1.0,qmin=2.0,V=[]):
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
    n=ncoord*(pts)
    mass_conv=mass*(amu/e_mass)
# In atomic units
    n1=n+1
    if dq==0.001:
        dq=(qmax-qmin)*((float(pts)+1.0)/(float(pts)-1.0))
    prefactor=(np.pi**2)/(4*mass_conv*dq**2)
    A=np.zeros((n,n),dtype=eval(numpy_precision))
# One has been added to i and j inside to make this consistent with paper
    for x in range (ncoord):
        for i in range(pts):
            for j in range(pts):
                if i==j:
                    A[i+x*pts,j+x*pts]=prefactor* (((2*n1**2+1)/3)-(np.sin(((i+1)*np.pi)/n1)**-2)) +V[i]
                else:
                    A[i+x*pts,j+x*pts]=prefactor* ((-1)**(i - j)) * ( np.sin((np.pi*(i-j)) / (2 * n1) )**-2  - 
                            np.sin((np.pi*(i+j+2)) / (2 * n1))**-2)
    #for x in A:
    #    print(x)
    return A
# I'll probably need to worry about the indicies this runs over when I attempt to generalize to multiple dimensions
#    dims=[]
#    for i in range(ncoord):
#        dims.append(pts)
#    indicies=np.zeros(dims)
# Slightly better indicies, in a list instead
#    from itertools import combinations_with_replacement
#    possiblevals=''
#    for i in range(pts):
#        possiblevals=possiblevals+str(i+1)
#    tmpindicies=combinations_with_replacement(possiblevals,ncoord)
#    indicies=[]
#    for x in tmpindicies:
#        tmpind2=[]
#        for y in range(len(x)):
#            tmpind2.append(int(x[y]))
#        indicies.append(tmpind2)
#    print(indicies)
# Idea to vectorize this function... not sure how to actually use position to evaluate, was thinking of setting up a sympy expression but gave up
#    A=np.identity(dimension,dtype=eval(numpy_precision))
#    A=np.piecewise(A,[A==0, A==1], [1, 1])
# See also http://stackoverflow.com/questions/21830112/evaluating-the-result-of-sympy-lambdify-on-a-numpy-mesgrid
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.mgrid.html

def main():
    constants(CODATA_year=2010)
    import numpy as np
    import sys
    if len(sys.argv)>1:
        potential=readpotential(sys.argv[1],r_units='angstrom')
    else:
        potential=readpotential(input('Give the file with the potential: '),r_units='angstrom')
    r_raw=np.array(potential[0],dtype=eval(numpy_precision))
    Energies_raw=np.array(potential[1],dtype=eval(numpy_precision))
    mass=potential[2]
    xmin,emin=returnsplinemin(r_raw,Energies_raw)
    r=r_raw-xmin
    Energies=Energies_raw-emin
    Ener_spline=cubicspline(r,Energies)
    xnew = np.linspace(min(r),max(r), num=num_points)
    vfit=returnsplinevalue(Ener_spline,xnew)
    Ham=H_array(ncoord=1,pts=len(r),mass=mass,V=Energies,qmax=max(r),qmin=min(r))
    eigenval, eigenvec=np.linalg.eig(Ham)
    Esort=np.sort(eigenval*hartreetocm)
# plotting stuff 
    if plotit:
        vfitcm=vfit*hartreetocm
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(r,Energies*hartreetocm,linestyle='none',marker='o')
        plt.plot(xnew,vfitcm,linestyle='solid',marker='None')
        mincut=[]
        maxcut=[]
        validE=[True]
        for x in range(1,len(Esort)-1,1):
            if np.multiply(0.95,(Esort[x]-Esort[x-1]))<(Esort[x+1]-Esort[x]):
                validE.append(True)
            else:
                Etoprint=x
                for y in range(x,len(Esort)+1):
                    validE.append(False)
                break
        for x in range(len(Esort)):
            if validE[x]:
                for y in range(len(xnew)):
                    if Esort[x]>vfitcm[y]:
                        mincut.append(xnew[y])
                        break
        for x in range(len(Esort)):
            if validE[x]:
                for y in range(len(xnew)-1,1,-1):
                    if Esort[x]>vfitcm[y]:
                        maxcut.append(xnew[y])
                        break
        for x in range(len(Esort)):
            if validE[x]:
                plt.plot((mincut[x],maxcut[x]),(Esort[x],Esort[x]),linestyle='solid')
#    plt.legend(['Points', 'Cubic Spline'])
        plt.title('Cubic-spline interpolation')
        plt.axis()
        plt.show()
    for x in range(Etoprint):
        print('{0:.{1}f}'.format(round(Esort[x+1]-Esort[x],num_print_digits),num_print_digits))


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

