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
        from numpy import multiply, divide, pi
        self.coordtype=int(input('coordinate type [r=0,theta=1,phi=2]: '))
        if self.coordtype==0:
            planckconstant=6.62606957E-34 # Js
            electron_charge=1.602176565E-19  #C
            rydberg= 109737.31568539 # cm-1
            light_speed= 299792458 # m/s
            bohr= (float(5) * light_speed * electron_charge**2)/ (planckconstant*rydberg) 
            self.maxval=float(input('Max value of r coordinate in bohr: '))
            self.minval=float(input('Min value of r coordinate in bohr: '))
            if r_unit==1:
                self.maxval=divide(self.maxval,bohr)
                self.minval=divide(self.minval,bohr)
        elif self.coordtype==1:
            self.maxval=pi
            self.minval=0.0
        elif self.coordtype==2:
            self.maxval=multiply(2,pi)
            self.minval=0.0
        self.mass=float(input('Give reduced mass: '))
        self.numpoints=int(input('number of points desired in coordinate: '))
        self.setgrid()

    def setgrid(self): 
        from numpy import linspace
        self.grid=linspace(self.minval,self.maxval,self.numpoints)

    def equivalencemwcoords(self,q_other):
        from scipy.optimize import minimize
        from numpy import subtract
        c=[]
        m1=self.mass
        m2=q_other.mass
        if self.coordtype==0 and q_other.coordtype==0:
            """ modify grid of two radial coordinates to make mass weighting equal"""
            c.append(self.maxval)
            c.append(self.minval)
            c.append(self.numpoints)
            c.append(q_other.maxval)
            c.append(q_other.minval)
            c.append(q_other.numpoints)
            bnds=((self.maxval*0.5,self.maxval*1.5),(self.minval*0.5,self.minval*1.5),(int(self.numpoints*0.1),self.numpoints*10),\
                    (q_other.maxval*0.5,q_other.maxval*1.5),(q_other.minval*0.5,q_other.minval*1.5),(int(q_other.numpoints*0.1),q_other.numpoints*10))
            result = minimize(frr, c, args=(m1,m2),bounds=bnds)
            if result.success:
                self.maxval, self.minval, self.numpoints, q_other.maxval, q_other.minval, q_other.numpoints = result.x
        elif self.coordtype==0:
            c.append(self.maxval)
            c.append(self.minval)
            c.append(self.numpoints)
            c.append(q_other.numpoints)
            bnds=((self.maxval*0.5,self.maxval*1.5),(self.minval*0.5,self.minval*1.5),(int(self.numpoints*0.1),self.numpoints*10),\
                    (int(q_other.numpoints*0.1),q_other.numpoints*10))
            q2range=subtract(q_other.maxval,q_other.minval)
            result = minimize(frang, c, args=(m1,m2,q2range),bounds=bnds)
            if result.success:
                self.maxval, self.minval, self.numpoints, q_other.numpoints = result.x
        elif q_other.coordtype==0:
            c.append(q_other.maxval)
            c.append(q_other.minval)
            c.append(q_other.numpoints)
            c.append(self.numpoints)
            bnds=((q_other.maxval*0.5,q_other.maxval*1.5),(q_other.minval*0.5,q_other.minval*1.5),(int(q_other.numpoints*0.1),q_other.numpoints*10),\
                    (int(self.numpoints*0.1),self.numpoints*10))
            q1range=subtract(self.maxval,self.minval)
            result = minimize(frang, c, args=(m1,m2,q1range),bounds=bnds)
            if result.success:
                q_other.maxval, q_other.minval, q_other.numpoints, self.numpoints  = result.x
        else:
            c.append(self.numpoints)
            c.append(q_other.numpoints)
            bnds=((int(self.numpoints*0.1),self.numpoints*10),(int(q_other.numpoints*0.1),q_other.numpoints*10))
            q2range=subtract(q_other.maxval,q_other.minval)
            q1range=subtract(self.maxval,self.minval)
            result = minimize(fangang, c, args=(m1,m2,q1range,q2range),bounds=bnds)
            if result.success:
                 self.numpoints, q_other.numpoints = result.x
        if result.success:
            print(result)
            self.setgrid()
            q_other.setgrid()
        else:
            print ('could not find a grid with equal mass weighting')

def con_pos(t):
    from numpy import subtract, abs
    return subtract(t,abs(t))

def frr(c,m1,m2):
    from numpy import subtract, power, round, divide, multiply
    mw1=multiply(power(divide(subtract(c[0],c[1]),subtract(round(c[2]),1)),2),m1)
    mw2=multiply(power(divide(subtract(c[3],c[4]),subtract(round(c[5]),1)),2),m2)
    return abs(subtract(mw2,mw1))

def frang(c,m1,m2,q2range):
    from numpy import subtract, power, round, divide, multiply
    mw1=multiply(power(divide(subtract(c[0],c[1]),subtract(round(c[2]),1)),2),m1)
    mw2=multiply(power(divide(q2range,subtract(round(c[3]),1)),2),m2)
    return abs(subtract(mw2,mw1))

def fangang(c,m1,m2,q1range,q2range):
    from numpy import subtract, power, round, divide, multiply
    mw1=multiply(power(divide(q1range,subtract(round(c[0]),1)),2),m1)
    mw2=multiply(power(divide(q2range,subtract(round(c[1]),1)),2),m2)
    return abs(subtract(mw2,mw1))

def massweightequal(dq1,m1,dq2,m2,printerrors=False):
    from numpy import subtract, power, multiply
    if abs(subtract(multiply(power(dq1,2),m1),multiply(power(dq2,2),m2)))> 1.0E-07:
        if printerrors:
            print('mass weighted spacing unequal as specified')
            print(dq1,m1,dq2,m2)
            print(m1*dq1**2,m2*dq2**2)
        return False
    return True

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
        if not massweightequal(dq1,coord[0].mass,dq2,coord[1].mass):
            coord[0].equivalencemwcoords(coord[1])
            dq1=np.subtract(coord[0].grid[1],coord[0].grid[0])
            dq2=np.subtract(coord[1].grid[1],coord[1].grid[0])
            if not massweightequal(dq1,coord[0].mass,dq2,coord[1].mass,printerrors=True):
                from sys import exit
                exit()
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
