#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 12:53:50 2016
pycharm
sublimetext
@author: zhencao
BITS = 10
"""
import math
import numpy as np

def SingleChain(M, N, K, L):
    xlo, xhi = -1.5, 33.5
    boxsize = 32
    ylo, yhi = 0.0, 32.0
    zlo, zhi = -2.0, 2.0 + K

    
    p = int((xhi - xlo) * 2)
    q = int(boxsize * (xhi -xlo) / 1.732) * 8
    r = int(boxsize * 5 / 1.732 * 8)
    s = int(boxsize * 2) #the number of beads in a line
    t = int(density * boxsize * K / chainlength)
    


    AtomTypes = 4
    BondTypes = 1
    AtomNo = L * M + 8 * q + 6 * r # +360*36+3600
    BondNo = (L - 1) * M  # +3600

    AtomList = []  # AtomId MolId AtomType charge x y z
    BondList = []  # BondIndex BondType A B

    for m in range(0, M):
        for n in range(1, N + 1):
            AtomList.append([n + m * L, m + 1, 1, m // t * 0.5 + 0.2 * (boxsize / chainlength) * n, m % t * 0.5, 0.15 * n])
            if (n != L):
                BondList.append([n + m * (chainlength - 1), 1, n + m * chainlength, n + 1 + m * chainlength])
        for n in range(N + 1, L + 1):
            AtomList.append([n + m * chainlength, m + 1, 2, m // t * 0.5 + 0.2 * (boxsize / chainlength) * n, m % t * 0.5, 0.15 * n])
            if (n != L):
                BondList.append([n + m * (chainlength - 1), 1, n + m * chainlength, n + 1 + m * chainlength])
    
    #set the upper and lower boundary
    for j in range(8):
        for i in range(j * q + 1, (j + 1) * q + 1):
            AtomList.append([M * chainlength + i, np.where(j < 4, M + 1, M + 2), np.where(j < 4, 3, 4), np.where(j % 2 == 0, i % p * 0.5 + (i // p) % 2 * 0.25, i % p * 0.5 + (i // p) % 2 * 0.25 + 0.25), np.where(j % 2 == 0, i // p * 0.43, i // p * 0.43 + 0.14), np.where(j < 4, -0.50 - 0.41 * j, K + 0.5 + (j - 4) * 0.41)])
    #set the side boundary
    for j in range(6):
        for i in range(j * r + 1, (j + 1) * r + 1):
            AtomList.append([M * chainlength + 8 * q + i, np.where(j < 3, M + 3, M + 4), 3, np.where(j < 3, -0.5 - j * 0.4, boxsize + 0.5 + (j - 3) * 0.4), np.where(j % 2 == 0, i % s * 0.5 + (i // s) % 2 * 0.25, i % s * 0.5 + (i // s) % 2 * 0.25 + 0.25), np.where(j % 2 == 0, (i % r) // s * 0.43, (i % r) // s * 0.43 + 0.14)])
    # post arrays
    #    for o in range(0, 20): #number of posts
    #      for k in range(0, int(r / 0.3)): #number of layers
    #               for j in range(1, 4):
    #                       for p in range(0, j*4):
    #                               if (j==1):
    #                                       AtomList.append([M*+33120*4+24*k+p+1+360*o,M+3,3,0.25*(j)*math.cos(2*3.14159/(j*4)*(p+(k%2)*0.5))+5.1+o%6*36.0,0.25*(j)*math.sin(2*3.14159/(j*4)*(p+(k%2)*0.5))+5.1+o//6*36.0,0.3*(k)])
    #                               if (j==2):
    #                                       AtomList.append([M*20+33120*4+24*k+4+p+1+360*o,M+3,3,0.25*(j)*math.cos(2*3.14159/(j*4)*(p+(k%2)*0.5))+5.1+o%6*36.0,0.25*(j)*math.sin(2*3.14159/(j*4)*(p+(k%2)*0.5))+5.1+o//6*36.0,0.3*(k)])
    #                               if (j==3):
    #                                       AtomList.append([M*20+33120*4+24*k+12+p+1+360*o,M+3,3,0.25*(j)*math.cos(2*3.14159/(j*4)*(p+(k%2)*0.5))+5.1+o%6*36.0,0.25*(j)*math.sin(2*3.14159/(j*4)*(p+(k%2)*0.5))+5.1+o//6*36.0,0.3*(k)])

    # brush layer
    #    for o in range(0, 36):
    #       for k in range(0, 5):
    #               for p in range(0, 4):
    #                       for j in range(0, 5):
    #                               AtomList.append([M*20+33120*4+360*36+100*o+20*k+j+p*5+1,M+4,2,(1.75+j)*math.cos(2*3.14159/(4)*(p+(k%2)*0.5))+5.1+o%6*36.0,(1.75+j)*math.sin(2*3.14159/(4)*(p+(k%2)*0.5))+5.1+o//6*36.0,0.9*(k)])
    #                               if (j == 0):
    #                                       BondList.append([19*M+100*o+20*k+j+(p*5)+1,1,M*20+33120*4+360*o+24*k*3+12+p*3+1,M*20+33120*4+360*36+100*o+20*k+j+p*5+1])
    #                               else:
    #                                       BondList.append([19*M+100*o+20*k+j+(p*5)+1,1,M*20+33120*4+360*36+100*o+20*k+j+p*5,M*20+33120*4+360*36+100*o+20*k+j+p*5+1])

    # writing headers and data into text file for simulation input
    f = open(str(M) + 'A' + str(N) + 'B' + str(chainlength - N) + '.dat', 'w')

    f.write('LAMMPS data file: Standing single chain in a box\n')
    f.write(' \n')
    f.write(str(AtomNo) + ' atoms\n')
    f.write(str(AtomTypes) + ' atom types\n')
    f.write(str(BondNo) + ' bonds\n')
    f.write(str(BondTypes) + ' bond types\n')
    f.write(' \n')#!/usr/bin/env python3
    f.write(str(xlo) + ' ' + str(xhi) + ' xlo xhi\n')
    f.write(str(ylo) + ' ' + str(yhi) + ' ylo yhi\n')
    f.write(str(zlo) + ' ' + str(zhi) + ' zlo zhi\n')
    f.write(' \n')
    f.write('Masses\n')
    f.write(' \n')
    for i in range(AtomTypes):
        f.write(str(i + 1) + ' 1\n')
    f.write(' \n')
    f.write('Atoms\n')
    f.write(' \n')
    for i in range(AtomNo):
        s = str(AtomList[i][0]) + " " + str(AtomList[i][1]) + " " + str(AtomList[i][2]) + " " + str(
            AtomList[i][3]) + " " + str(AtomList[i][4]) + " " + str(AtomList[i][5])
        s = s + '\n'
        f.write(s)

    f.write('\nBonds\n')
    f.write(' \n')

    for i in range(BondNo):
        s = str(BondList[i][0]) + " " + str(BondList[i][1]) + " " + str(BondList[i][2]) + " " + str(BondList[i][3])
        s = s + '\n'
        f.write(s)

    f.close()


N = [2]#number of A bead
boxsize = 32
chainlength = 10
density = 5
k = 5
for n in N:
    SingleChain(int(k * boxsize**2 * density / chainlength), n, k, chainlength) #43 = box size ** 2 / chain length

