#!/usr/bin/env python3
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (bytes, str, open, super, range, zip, round, input, int, pow, object)
class q:
    r_unit=-1

    def __init__(self):
        r_unit=0
#        try:
#            r_unit
#        except NameError:
#            r_unit=int(input('Units 0=a.u. 1=angstrom: '))
        import numpy as np
        self.coordtype=int(input('coordinate type [r=0,theta=1,phi=2]: '))
        if self.coordtype==0:
            planckconstant=6.62606957E-34 # Js
            electron_charge=1.602176565E-19  #C
            rydberg= 109737.31568539 # cm-1
            light_speed= 299792458 # m/s
            bohr= (float(5) * light_speed * electron_charge**2)/ (planckconstant*rydberg) 
            self.maxval=float(input('Max value of r coordinate: '))
            self.minval=float(input('Min value of r coordinate: '))
            if r_unit==1:
                self.maxval=np.divide(self.maxval,bohr)
                self.minval=np.divide(self.minval,bohr)
        elif self.coordtype==1:
            self.maxval=np.pi
            self.minval=0.0
        elif self.coordtype==2:
            self.maxval=np.multiply(2,np.pi)
            self.minval=0.0
        self.mass=float(input('Give reduced mass: '))
        self.numpoints=int(input('number of points desired in coordinate: '))
        self.grid=np.linspace(self.minval,self.maxval,self.numpoints)

def massweightequal(dq1,m1,dq2,m2):
    import numpy as np
    if np.abs(np.subtract(np.multiply(np.power(dq1,2),m1),np.multiply(np.power(dq2,2),m2)))> 1.0E-07:
        from sys import exit
        exit('mass weighted grid spacings are unequal')

def main():
    fieldlength=18
    outfile='tmp.pot'
    numcoordinates=int(input('number of coordinates: '))
    coord=[]
    for x in range(int(numcoordinates)):
        coord.append(q())
    import numpy as np
    coordtypedict={0: 'radial',1: 'angular_pi', 2: 'angular_2pi'}
    gridout=[]
    if numcoordinates==1:
        indices=np.array(np.meshgrid(np.arange(coord[0].numpoints))).T.reshape(-1,1)
        for x in range(len(indices)):
            gridout.append(str(coord[0].grid[indices[x][0]]).ljust(fieldlength))
    elif numcoordinates==2:
        dq1=np.subtract(coord[0].grid[1],coord[0].grid[0])
        dq2=np.subtract(coord[1].grid[1],coord[1].grid[0])
        massweightequal(dq1,coord[0].mass,dq2,coord[1].mass)
        indices=np.array(np.meshgrid(np.arange(coord[0].numpoints),np.arange(coord[1].numpoints))).T.reshape(-1,2)
        indices=np.split(indices,numcoordinates,axis=1) 
        for x in range(len(indices[0])):
            gridout.append(str(coord[0].grid[indices[0][x][0]]).ljust(fieldlength) + str(coord[1].grid[indices[1][x][0]]).ljust(fieldlength))
    from os.path import isfile
    if isfile(outfile):
        if 'y' not in input('outfile exists, overwrite? [y,N]').lower():
            from sys import exit
            exit()
    txt=open(outfile,'w')
    txt.write('mass= ')
    for x in range(numcoordinates):
        txt.write(str(coord[x].mass).ljust(fieldlength) + ' ')
    txt.write('\n')
    if q.r_unit==1:
        txt.write('angstrom\n')
    else:
        txt.write('bohr\n')
    for x in range(numcoordinates):
        txt.write(str(coordtypedict[coord[x].coordtype]) + ' ')
    for x in range(len(gridout)):
        txt.write('\n'+gridout[x])

if __name__=="__main__":
    main()
