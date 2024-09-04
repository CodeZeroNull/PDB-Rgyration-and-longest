# -*- coding: utf-8 -*-
#!/usr/bin/python

import math 
import sys
import numpy as np
from scipy.spatial.distance import pdist, squareform


def RgLong(filename):
    '''
    Given a PDB structure file, calculates the Radius of Gyration (Rg) and
    the longest distance.
    NB: Only takes in account H, C, O, N and S atoms.
    
    Returns (values in Angstrom):
        - the Rg
        - coordinates for the center of mass
        - the longest distance between atoms
        - coordinates for the first atom furthest apart
        - coordinates for the second atom furthest apart
    '''


    # Get coordinates and masses for each atom,
    # ignoring lines that do not end with: H, C, O, N, or S
    coord = list()
    mass = list()
    Structure = open(filename, 'r')
    for line in Structure:
        try:
            line = line.split()
            # Owing to split failing to get correct coordinates for some
            # alternate conformations where space is removed, I am going to
            # make a quick and easy and nasty hack here, ideally change things
            # so that files are imported using the pdbreader library instead
            if len(line) == 13:
                x = float(line[6])
                y = float(line[7])
                z = float(line[8])
                occ = float(line[9])
            elif len(line) == 12:
                x = float(line[5])
                y = float(line[6])
                z = float(line[7])
                occ = float(line[8])
            coord.append([x, y, z])
            if line[-1] == 'C':
                mass.append(12.0107*occ)
            elif line[-1] == 'O':
                mass.append(15.9994*occ)
            elif line[-1] == 'O1-':
                mass.append(15.9994*occ)
            elif line[-1] == 'N':
                mass.append(14.0067*occ)
            elif line[-1] == 'N1+':
                mass.append(14.0067*occ)
            elif line[-1] == 'S':
                mass.append(32.065*occ)
            elif line[-1] == 'H':
                mass.append(1.00794*occ)
        except:
            pass
    Structure.close()
    
    
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
    CM = [cmx, cmy, cmz]
        
    # Calculate Radius of Gyration
    rg = 0
    for (i, j, k), m in zip(coord, mass):
        rg+=m*((i-cmx)**2+(j-cmy)**2+(k-cmz)**2)
    rg/=tmass
    rg=math.sqrt(rg)


    # Determine the Longest distance in the structure
    
    # Convert distance vector to a condensed distance matrix, then
    # converts that matrix  into a redundant distance matrix
    # (size is 'len(coord) x len(coord)')
    distances = squareform(pdist(coord))
    
    # Find maxmimum distance and index of maximum distance
    Longest = np.max(distances) 
    indexLongest = np.argmax(distances)
    
    #Find coordinates of atoms at longest distance
    indexLrow, indexLcol = np.unravel_index(indexLongest, distances.shape)
  
     
    return(round(rg, 1), CM, round(Longest, 1), coord[indexLrow], coord[indexLcol])



if __name__ == '__main__':
    Rg,CM,Long,point1,point2 = RgLong(sys.argv[1])
    print("Radius of Gyration is ", Rg)
    print("Longest distance is ", Long)
    print()
    print("One can check the radius of gyration in PyMOL with the command:")
    print("pseudoatom tmpPointCM, pos=", CM, ", vdw=", Rg)
    print("Remeber to select tmpPoint and show as Spheres and maybe 'set sphere_transparency=0.5' ")
    print()
    print("One can get the distance in PyMOL with the commands:")
    print("pseudoatom tmpPoint1, pos=", point1)
    print("pseudoatom tmpPoint2, pos=", point2)
    print("dist distance1, tmpPoint1, tmpPoint2")
    print()
