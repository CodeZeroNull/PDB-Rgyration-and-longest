# -*- coding: utf-8 -*-
#!/usr/bin/python

import math 
import sys

def Rg(filename):
    '''
    Given a PDB structure file, calculates the Radius of Gyration (Rg) and
    the longest distance.
    NB: Only takes in account H, C, O, N and S atoms.
    
    Returns the Rg and Longest values in Angstrom.
    '''
    # Get coordinates and masses for each atom,
    # ignoring lines that do not end with: H, C, O, N, or S
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
    
    # Determine Center of Mass Coordinates and Total Mass
    cmx,cmy,cmz,tmass=(0.0,0.0,0.0,0.0)
    for (i, j, k), m in zip(coord, mass):
        tmass+=m
        cmx+=m*i
        cmy+=m*j
        cmz+=m*k
    cmx/=tmass
    cmy/=tmass
    cmz/=tmass
    print("One can check the position in PyMOL with the command:")
    print("pseudoatom tmpPoint, pos=[", cmx, "," , cmy, "," , cmz, "]")
        
    # Calculate Radius of Gyration
    rg = 0
    for (i, j, k), m in zip(coord, mass):
        rg+=m*((i-cmx)**2+(j-cmy)**2+(k-cmz)**2)
    rg/=tmass
    rg=math.sqrt(rg)
    print("One can check the radius of gyration in PyMOL with the command:")
    print("pseudoatom tmpPoint, pos=[", cmx, "," , cmy, "," , cmz, "], vdw=", rg)
    print("Remeber to select tmpPoint and show as Spheres")

    # Determine Longest distance in the structure
    Longest = 0
     
    return(round(rg, 1 ), Longest)

if __name__ == '__main__':
    print('Rg = {}'.format(Rg(sys.argv[1])))
