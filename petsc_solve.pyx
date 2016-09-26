cimport numpy as np
def H_array_petsc(pts=5,coordtype=['r'],mass=[0.5],qmin=[1.0],qmax=[2.0],V=[]):
    import sys, slepc4py
    slepc4py.init(sys.argv)
    from petsc4py import PETSc
    from slepc4py import SLEPc
    cdef int totpts=1
    cdef int i, j, k, iter
    cdef int [:,:,:,:] indices
    cdef int [:,:] ijindex
    for x in range(len(pts)):
        totpts=totpts*pts[x]
    opts=PETSc.Options()
    from mpi4py import MPI
    import numpy as np
    comm = MPI.COMM_WORLD
    cdef int rank = comm.Get_rank()
    A=PETSc.Mat().create()
    A.setSizes(totpts,totpts)
    A.setFromOptions()
    A.setUp()
    rstart,rend=A.getOwnershipRange()
    """  Note the ownership seems to entirely be localted on one process..."""
    ncoord=len(coordtype)
    if rank==0:
        """ Following https://pythonhosted.org/slepc4py/usrman/tutorial.html
        note that on a 10^4 matrix the construction takes 22 minutes, solving for first 3 eigenvalues is 1min
        solving for all is 86 minutes. For optimizer I'll need to improve assembly time as well..."""
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
#            indices=np.array(np.meshgrid(np.arange(ptspercoord),np.arange(ptspercoord),np.arange(ptspercoord),np.arange(ptspercoord)\
#                    ,np.arange(ptspercoord),np.arange(ptspercoord))).T.reshape(-1,6)
        indices=np.array(np.split(indices,ncoord*2,axis=1),dtype=eval(inttype))
        k=0
        if ncoord==1:
            for i in range(totpts):
                for j in range(totpts):
                    if i==j:
                        A[i,i] =np.add(D1[i,j],V[i])
                    else:
                        A[i,j]=D1[i,j]
        elif ncoord==2:
            totiter=(totpts)**2
            ijindex=np.zeros((totiter,2),dtype=np.uint64)
            iter=0
            for i in range(totpts):
                for j in range(totpts):
                    ijindex[iter,0]=i
                    ijindex[iter,1]=j
                    iter+=1
            iter=0
            while iter<totiter:
                if np.equal(indices[0,iter],indices[3,iter]):
                    if np.equal(indices[1,iter],indices[2,iter] ):
                        A[ijindex[iter]]=np.add(np.add(D1[indices[0,iter],indices[3,iter]], D2[indices[1,iter],indices[2,iter]]),V[k])
                        k+=1
                    else:
                        A[ijindex[iter]]=D2[indices[1,iter],indices[2,iter]]
                elif np.equal(indices[1,iter],indices[2,iter] ): 
                    A[ijindex[iter]]=D1[indices[0,iter],indices[3,iter]]
                    for i in range(pts[1]-1-np.squeeze(indices[1,iter])):
                        iter+=1
                iter+=1
        else:
            raise ValueError("not implemented for %d dimensions" % (ncoord))
    A.assemble()
    E = SLEPc.EPS(); E.create()
    E.setOperators(A)
    E.setProblemType(SLEPc.EPS.ProblemType.HEP)
    E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)
#    E.setDimensions(6)
    E.setDimensions(totpts)
    E.setFromOptions()
    E.solve()
#    E.view()



    Print = PETSc.Sys.Print

    its = E.getIterationNumber()
#    Print("Number of iterations of the method: %d" % its)
    eps_type = E.getType()
#    Print("Solution method: %s" % eps_type)
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

def H_array_1d(pts=5,coordtype='r',mass=0.5,qmin=1.0,qmax=2.0):
    import numpy as np
    n=pts
    numpy_precision="np.float64"
    A=np.zeros((n,n),dtype=eval(numpy_precision))
# In atomic units
# One has been added to i and j inside to make this consistent with paper
    amu=1.660538921e-27
    e_mass= 9.10938215e-31
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

def mwspace(coordtype='r',rlen=1.0,mass=1.0,pts=2):
    from numpy import multiply, power, divide, pi , add
    if coordtype=='r':
        mwspace=multiply(power(divide(rlen,add(pts,1)),2),mass)
    elif coordtype=='theta':
        mwspace=multiply(power(divide(pi,add(pts,1)),2),mass)
    elif coordtype=='phi':
        mwspace=multiply(power(divide(multiply(2,pi),pts),2),mass)
    return mwspace

#cdef int num_print_digits=3
#cdef double light_speed= 299792458
#cdef double planckconstant=6.62606957e-34
#cdef double amu=1.660538921e-27
#cdef double e_mass= 9.10938215e-31
#cdef double electron_charge=1.602176565e-19
#cdef double rydberg= 109737.31568539
#cdef double bohr= (float(5) * light_speed * electron_charge**2)/ (planckconstant*rydberg) 
#cdef double hartreetocm= 2*rydberg


