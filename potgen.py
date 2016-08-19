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
        from numpy import linspace, delete
        self.grid=linspace(self.minval,self.maxval,self.numpoints)
        #if self.coordtype==0: 
            #Actual grid here doesn't include a or b for r 
            #Note for gridspacing this has no effect, it simply doesn't include the points.  
            #Reason for this is the calculation is the potential is not included at a or b 
        if  self.coordtype==1:
            #Actual grid here doesn't include 0 or pi for theta"""
            self.grid=delete(self.grid,0)
            self.grid=delete(self.grid,len(self.grid)-1)
            self.numpoints=self.numpoints-2
        elif self.coordtype==2:
            #Actual grid here doesn't include 2pi for phi
            self.grid=linspace(self.minval,self.maxval,self.numpoints)
            self.grid=delete(self.grid,len(self.grid)-1)
            self.numpoints=self.numpoints-1

    def equivalencemwcoords(self,q_other,forcegrid=True):
        from scipy.optimize import minimize, basinhopping
        from numpy import subtract, multiply, int, pi, add
        m1=self.mass
        m2=q_other.mass
        trialval=[]
        triali=[]
        trialj=[]
        vals=[]
        if self.coordtype==0 and q_other.coordtype==0 and not forcegrid:
            """ modify grid of two radial coordinates to make mass weighting equal"""
            c=[self.maxval, self.minval, q_other.maxval, q_other.minval]
            minbnds=(multiply(self.maxval,0.9),multiply(self.minval,0.9),\
                    multiply(q_other.maxval,0.9),multiply(q_other.minval,0.9))
            maxbnds=(multiply(self.maxval,1.1),multiply(self.minval,1.1),\
                multiply(q_other.maxval,1.1),multiply(q_other.minval,1.1))
            bnds=list(zip(minbnds,maxbnds))
            for i in range(int(round(self.numpoints*0.5)),5*self.numpoints, 1):
                for j in range(int(round(q_other.numpoints*0.5)),5*q_other.numpoints, 1):
                    result=minimize(frr,c,args=(m1,m2,i,j),bounds=bnds)
                    if result.fun<1e-07:
                        trialval.append(result.fun)
                        triali.append(i)
                        trialj.append(j)
                        vals.append(result.x)
        elif self.coordtype==0 and not forcegrid:
            """ modify grid of one radial coordinates to make mass weighting equal"""
            c=[self.maxval, self.minval]
            q2range=subtract(q_other.maxval,q_other.minval)
            q1range=subtract(self.maxval,self.minval)
            minbnds=(multiply(self.maxval,0.9),multiply(self.minval,0.9))
            maxbnds=(multiply(self.maxval,1.1),multiply(self.minval,1.1))
            bnds=list(zip(minbnds,maxbnds))
            for i in range(int(round(self.numpoints*0.5)),5*self.numpoints, 1):
                for j in range(int(round(q_other.numpoints*0.5)),5*q_other.numpoints, 1):
                    result=minimize(frang,c,args=(m1,m2,q2range,i,j),bounds=bnds)
                    if result.fun<1e-07:
                        trialval.append(result.fun)
                        triali.append(i)
                        trialj.append(j)
                        vals.append(result.x)
        elif q_other.coordtype==0 and not forcegrid:
            """ modify grid of one radial coordinates to make mass weighting equal"""
            c=[q_other.maxval ,q_other.minval]
            minbnds=(multiply(q_other.maxval,0.9),multiply(q_other.minval,0.9))
            maxbnds=(multiply(q_other.maxval,1.1),multiply(q_other.minval,1.1))
            q1range=subtract(self.maxval,self.minval)
            bnds=list(zip(minbnds,maxbnds))
            for i in range(int(round(self.numpoints*0.5)),5*self.numpoints, 1):
                for j in range(int(round(q_other.numpoints*0.5)),5*q_other.numpoints, 1):
                    """ i and j are reversed to keep things consistent below"""
                    result=minimize(frang,c,args=(m2,m1,q1range,j,i),bounds=bnds)
                    if result.fun<1e-07:
                        trialval.append(result.fun)
                        triali.append(i)
                        trialj.append(j)
                        vals.append(result.x)
        else:
            """ modify number of points in grid to make mass weighting equal"""
            q1range=subtract(self.maxval,self.minval)
            q2range=subtract(q_other.maxval,q_other.minval)
            for i in range(int(round(self.numpoints*0.3)),10*self.numpoints, 1):
                for j in range(int(round(q_other.numpoints*0.3)),10*q_other.numpoints, 1):
                    val1= fangang(m1,m2,q1range,q2range,i,j)
                    if val1<1e-07:
                        trialval.append(val1)
                        triali.append(i)
                        trialj.append(j)
        if len(triali)>1:
            if self.coordtype==0 and q_other.coordtype==0 and not forcegrid:
                print('{0:5} {1:4} {2:4} {3:15} {4:15}'\
                        .format('Number','pts1','pts2','q1','q2'))
                for i in range(len(triali)):
                    print('{0:5} {1:5} {2:5} {4:.3f}-{3:.3f} {6:.3f}-{5:.3f}'\
                            .format(i,triali[i],trialj[i],vals[i][0],vals[i][1],vals[i][2],vals[i][3]))
                itouse=int(input('choose the number corresponding to the desired number of points: '))
                self.maxval, self.minval, q_other.maxval, q_other.minval  = vals[itouse]
            elif self.coordtype==0 and not forcegrid:
                print('{0:5} {1:4} {2:4} {3:15} '.format('Number','pts1','pts2','q1'))
                for i in range(len(triali)):
                    print('{0:5} {1:5} {2:5} {4:.3f}-{3:.3f}'.format(i,triali[i],trialj[i],vals[i][0],vals[i][1]))
                itouse=int(input('choose the number corresponding to the desired number of points: '))
                self.maxval, self.minval  = vals[itouse]
            elif q_other.coordtype==0 and not forcegrid:
                print('{0:5} {1:4} {2:4} {3:15} '.format('Number','pts1','pts2','q2'))
                for i in range(len(triali)):
                    print('{0:5} {1:5} {2:5} {4:.3f}-{3:.3f}'.format(i,triali[i],trialj[i],vals[i][0],vals[i][1]))
                itouse=int(input('choose the number corresponding to the desired number of points: '))
                q_other.maxval, q_other.minval = vals[itouse]
            else:
                print('{0:5} {1:4} {2:4}'.format('Number','pts1','pts2'))
                for i in range(len(triali)):
                    print('{0:5} {1:5} {2:5}'.format(i,triali[i],trialj[i]))
                itouse=int(input('choose the number corresponding to the desired number of points: '))
            self.numpoints=triali[itouse]
            q_other.numpoints=trialj[itouse]
        else:
            from sys import exit
            print('no reasonable solutions found for this potential')
            exit()
        self.setgrid()
        q_other.setgrid()

def frr(c,m1,m2,pts1,pts2):
    from numpy import subtract, power, round, divide, multiply#, add, int, abs
    mw1=multiply(power(divide(subtract(c[0],c[1]),subtract((pts1),1)),2),m1)
    mw2=multiply(power(divide(subtract(c[2],c[3]),subtract((pts2),1)),2),m2)
    return abs(subtract(mw2,mw1))
#    return add(abs(subtract(mw2,mw1)),multiply(0.01,add(abs(subtract(int(c[2]),c[2])),abs(subtract(int(c[5]),c[5])))))

def frang(c,m1,m2,q2range,pts1,pts2):
    from numpy import subtract, power, round, divide, multiply
    mw1=multiply(power(divide(subtract(c[0],c[1]),subtract((pts1),1)),2),m1)
    mw2=multiply(power(divide(q2range,subtract((pts2),1)),2),m2)
    return abs(subtract(mw2,mw1))

def fangang(m1,m2,q1range,q2range,pts1,pts2):
    from numpy import subtract, power, round, divide, multiply
    mw1=multiply(power(divide(q1range,subtract(pts1,1)),2),m1)
    mw2=multiply(power(divide(q2range,subtract(pts2,1)),2),m2)
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
    modifiedinputfile='tmp.inp'
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
        if 'y' not in input('outfile (tmp.pot) exists, overwrite? [y,N]').lower():
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
    txt.write('\n')
    txt.write('---\n')
    for x in range(numcoordinates):
        coordset=sorted(set(coord[x].grid))
        outstr='q'+str(x)+'(1:phh)=['
        for y in range(len(coordset)): 
            if len(outstr)<930 and y!=(len(coordset)-1):
                outstr=outstr+str(coordset[y])+', '
            elif len(outstr)>930:
                txt.write(outstr[:-2].replace('phh',str(y))+']\n')
                outstr='q'+str(numcoordinates)+'('+str(y+1)+':phh)=['
                outstr=outstr+str(coordset[y])+', '
            elif y==(len(coordset)-1):
                outstr=outstr+str(coordset[y])
                txt.write(outstr.replace('phh',str(y))+']\n')
        txt.write('\n')
    txt.close()
    if isfile(modifiedinputfile):
        if 'y' not in input('input file (tmp.inp) exists, overwrite? [y,N]').lower():
            from sys import exit
            exit()
    txt=open(modifiedinputfile,'w')
    txt.write(str(numcoordinates)+'\n')
    for x in range(numcoordinates):
        txt.write(str(coord[x].coordtype)+'\n')
        if coord[x].coordtype==0:
            txt.write(str(coord[x].maxval)+'\n')
            txt.write(str(coord[x].minval)+'\n')
        txt.write(str(coord[x].mass)+'\n')
        if coord[x].coordtype==0:
            txt.write(str(coord[x].numpoints)+'\n')
        elif coord[x].coordtype==1:
            txt.write(str(coord[x].numpoints+2)+'\n')
        elif coord[x].coordtype==2:
            txt.write(str(coord[x].numpoints+1)+'\n')

if __name__=="__main__":
    main()

