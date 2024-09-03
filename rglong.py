# -*- coding: utf-8 -*-
#!/usr/bin/python

import math 
import sys

def Rg(filename):
    '''
    Calculates the Radius of Gyration (Rg) of a protein given its .pdb 
    structure file.
    NB: Only takes in account H, C, O, N and S atoms.
    
    Returns the Rg value in Angstrom.
    '''
    # Get coordinates and masses for each atom,
    # ignoring lines that do not end with either: H, C, O, N, S
    coord = list()
    mass = list()
    Structure = open(filename, 'r')
    for line in Structure:
        try:
            line = line.split()
            x = float(line[6])
            y = float(line[7])
            z = float(line[8])
            coord.append([x, y, z])
            if line[-1] == 'C':
                mass.append(12.0107)
            elif line[-1] == 'O':
                mass.append(15.9994)
            elif line[-1] == 'N':
                mass.append(14.0067)
            elif line[-1] == 'S':
                mass.append(32.065)
            elif line[-1] == 'H':
                mass.append(1.00794)
        except:
            pass
    
    # Determine Center of Mass
    
    # Calculate Radius of Gyration
    xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
    tmass = sum(mass)
    rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
    mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
    rg = math.sqrt(rr / tmass-mm)
    return(round(rg, 1 ))

if __name__ == '__main__':
    print('Rg = {}'.format(Rg(sys.argv[1])))
