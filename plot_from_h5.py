def ploteigenfunctions():
    constants(CODATA_year=2010)
    import numpy as np
    import sys
    if len(sys.argv)>1 and 'h5' in sys.argv[1]:
        from dvr import loadeigen
        if len(sys.argv)>2:
            loadeigen(eigfile=sys.argv[1],eigenvectoplot=int(sys.argv[2]))
        else:
            loadeigen(eigfile=sys.argv[1],eigenvectoplot=int(input('number of eigenvectors to plot:'))) 
    else:
        print('h5 input file not specified properly to print eigenfunction')

if __name__=="__main__":
    ploteigenfunctions()
