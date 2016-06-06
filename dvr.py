#!/usr/bin/env python3
""" This program was written following 
Citation: The Journal of Chemical Physics 96, 1982 (1992); doi: 10.1063/1.462100
as a guide
Should use  scipy.sparse.linalg.eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')[source] as a way of doing Lanczos...."""
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (bytes, str, open, super, range, zip, round, input, int, pow, object)

def constants(CODATA_year=2010):
    """ CODATA constants used in program and other global defintions including number of points and numpy datatype for reals"""
    global numpy_precision, num_points
    numpy_precision="np.float64"
    num_points=1000
    global planckconstant, light_speed, Rydberg, electron_charge, amu, bohr
    light_speed= 299792458 # m/s
    if CODATA_year==2010:
        planckconstant=6.62606957E-34 # Js
        amu=1.660538921E-27 #kg
        electron_charge=1.602176565E-19  #C
        rydberg= 109737.31568539 # cm-1
    elif CODATA_year==2006:
        planckconstant= 6.626068963E-34 # Js
        amu= 1.660538782E-27 # kg
        electron_charge= 1.602176487E-19 # C
        rydberg= 109737.31568527 #cm-1
    else:
        from sys import exit
        exit('Constants not found')
    bohr= (float(5) * light_speed * electron_charge**2)/ (planckconstant*rydberg) 
    hartree= 2*rydberg


def openandread(filename):
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
    # potential should have coordinate and units as main input
    # anyline starting with ! or # is commented out
    commentoutchars=['!','#']
    lines=openandread(inp)
    r=[]
    energy=[]
    if r_units.lower()=='bohr':
        r_unitconversion=1.0
    elif r_units.lower()=='angstorm':
        r_unitconversion=(1.0/bohr)
#        r_unitconversion=1.88972613
    else:
        from sys import exit
        exit('No valid units given for length')
    for x in lines:
        if x[0] not in commentoutchars:
            linesplit=x.split()
            if len(linesplit)==2:
                r.append(float(linesplit[0])*r_unitconversion)
                energy.append(float(linesplit[1]))
    return (r,energy)


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
    from scipy import interpolate
    return interpolate.splrep(x, y, s=0)

def returnsplinevalue(spline,xnew):
    """ returns the value(s) of a spline function at a given x(s)"""
    from scipy import interpolate
    return interpolate.splev(xnew, spline, der=0)

def T_array(n=1,reduced_mass=0.5,dx=0.001):
    """ Kinetc Energy Array (dimensionality=2): see Eq A6a and A6b of JCP 96, 1982 (1992): note 
    constants are defined in constants module globally earlier"""
    import numpy as np
#    A=np.identity(dimension,dtype=eval(numpy_precision))
    prefactor=(np.pi**2*(planckconstant/(2*np.pi))**2/(4*reduced_mass*(dx**2)))
    A=np.zeros((n,n),dtype=eval(numpy_precision))
    for i in range(n):
        for j in range(n):
            if i==j:
                A[i,j]=prefactor*(((2.0*num_points**2+1.)/3.0)-(1./np.sin((np.pi*i)/num_points)))
            else:
                A[i,j]=prefactor*
    return A

#    A=np.piecewise(A,[A==0, A==1], [1, 1])

def main():
    constants(CODATA_year=2010)
    import numpy as np
    import sys
    T=T_array(n=5,reduced_mass=(1./3.),dx=0.001)
    from sys import exit
    print(T)
    exit()
    if len(sys.argv)>1:
        potential=readpotential(sys.argv[1])
    else:
        potential=readpotential(input('Give the file with the potential: '))
    r=np.array(potential[0],dtype=eval(numpy_precision))
    Energies=np.array(potential[1],dtype=eval(numpy_precision))
    Ener_spline=cubicspline(r,Energies)
    xnew = np.linspace(min(r),max(r), num=num_points)
    ynew=returnsplinevalue(Ener_spline,xnew)
# plotting stuff to make sure splines work
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(r,Energies,linestyle='none',marker='o')
    plt.plot(xnew,ynew,linestyle='solid',marker='None')
#    plt.legend(['Points', 'Cubic Spline'])
    plt.title('Cubic-spline interpolation')
    plt.axis()
    plt.show()


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

