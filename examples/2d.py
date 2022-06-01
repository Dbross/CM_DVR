#!/usr/bin/env python3
def main():
    import sys
    if len(sys.argv)>1:
        txt=open(sys.argv[1])
    else:
        txt=open(input('Give the input potential to double'))
    lines=txt.readlines()
    txt.close()
    E=[]
    r=[]
    for x in lines:
        splitline=x.split()
        r.append(splitline[0])
        E.append(splitline[1])
    r1=[]
    r2=[]
    Efin=[]
    for x in range(len(r)):
        for y in range(len(r)):
            r1.append(r[x])
            r2.append(r[y])
            Efin.append(float(E[x])+float(E[y]))
    for z in range(len(r1)):
        print (r1[z],r2[z],Efin[z])



if __name__=="__main__":
    main()

