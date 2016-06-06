#!/usr/bin/env python3
# This program was written following 
# Citation: The Journal of Chemical Physics 96, 1982 (1992); doi: 10.1063/1.462100
# as a guide
# Should use  scipy.sparse.linalg.eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')[source] as a way of doing Lanczos....

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
    if r_units.lower=='bohr':
        r_unitconversion=1.0
    elif r_units.lower=='angstorm':
        r_unitconversion=1.88972613
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

def main():
    from numpy import array
    A = array([[2.0,1.0],[5.0,7.0]])
    b = array([11.0,13.0])
    guess = array([1.0,1.0])

    sol = jacobi(A,b,N=25,x=guess)
    from pprint import pprint
    pprint(sol)
    from sys import exit
    exit()
    import sys
    if len(sys.argv)>1:
        potential=readpotential(sys.argv[1])
    else:
        potential=readpotential(input('Give the file with the potential: '))


if __name__=="__main__":
    main()

