
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


def SingleChain(M, N, K, box, chain_length, post_array):
    xlo, xhi = box[0], box[1]
    ylo, yhi = box[2], box[3]
    zlo, zhi = box[4], box[5]
    ysize = yhi - ylo
    post_number = post_array.shape[0]
    radium = 0.15
    bond2_length = 0.2
    #p = int((xhi - xlo) * 2)
    q = int(xsize * (xhi - xlo) / 1.732) * 8
    #r = int(xsize * K * 0.69 / 1.732 * 8)
    #s = int(xsize * 2)  # the number of beads in a line
    t = int(K / 0.3)  # the number of the layer

    AtomTypes = 3
    BondTypes = 2
    Post_array_number = post_number * t * 12
    AtomNo = chain_length * M + Post_array_number + 4 * t * post_number  # 6 * r # +360*36+3600
    BondNo = (chain_length - 1) * M + 4 * t * post_number#+3600
    AtomList = []  # AtomId MolId AtomType charge x y z
    BondList = []  # BondIndex BondType A B

    # lay particles
    a = int(density * ysize * z / chain_length * 2)
    b = 0.9
    for m in range(0, M):
        for n in range(1, N + 1):

            AtomList.append([n + m * chain_length, m + 1, 1, m // a + b * n, m % a, 0.15 * n])
            if (n != chain_length):
                BondList.append([n + m * (chain_length - 1), 1, n + m * chain_length, n + 1 + m * chain_length])
        for n in range(N + 1, chain_length + 1):
            AtomList.append([n + m * chain_length, m + 1, 2, m // a + b * n, m % a, 0.15 * n])
            if (n != chain_length):
                BondList.append([n + m * (chain_length - 1), 1, n + m * chain_length, n + 1 + m * chain_length])

    # set the upper and lower boundary
    #for j in range(8):
    #    for i in range(j * q + 1, (j + 1) * q + 1):
    #        AtomList.append([M * chain_length + i, np.where(j < 4, M + 1, M + 2), np.where(j < 4, 3, 4),
    #                         np.where(j % 2 == 0, i % p * 0.5 + (i // p) % 2 * 0.25,
    #                                  i % p * 0.5 + (i // p) % 2 * 0.25 + 0.25),
    #                         np.where(j % 2 == 0, i // p * 0.43, i // p * 0.43 + 0.14),
    #                         np.where(j < 4, -0.50 - 0.41 * j, K + 0.5 + (j - 4) * 0.41)])
    # set the side boundary
    # for j in range(6):
    #    for i in range(j * r + 1, (j + 1) * r + 1):
    #        AtomList.append([M * 4 + 8 * q + i, np.where(j < 3, M + 3, M + 4), 3, np.where(j < 3, -0.5 - j * 0.4, xsize + 0.5 + (j - 3) * 0.4), np.where(j % 2 == 0, i % s * 0.5 + (i // s) % 2 * 0.25, i % s * 0.5 + (i // s) % 2 * 0.25 + 0.25), np.where(j % 2 == 0, (i % r) // s * 0.43, (i % r) // s * 0.43 + 0.14)])
    # post arrays
    print(post_number)
    for o in range(post_number):  # number of posts
        for k in range(t):  # number of layers
            for j in range(1, 3):
                for p in range(0, j * 4):
                    if (j == 1):
                        AtomList.append([M * chain_length + 12 * k + p + 1 + 12 * t * o, M + 2, 3,
                                         radium * (j) * math.cos(2 * 3.14159 / (j * 4) * (p + (k % 2) * 0.5)) +
                                         post_array[o][0],
                                         radium * (j) * math.sin(2 * 3.14159 / (j * 4) * (p + (k % 2) * 0.5)) +
                                         post_array[o][1], 0.3 * (k)])
                    if (j == 2):
                        AtomList.append([M * chain_length + 12 * k + 4 + p + 1 + 12 * t * o, M + 2, 3,
                                         radium * (j) * math.cos(2 * 3.14159 / (j * 4) * (p + (k % 2) * 0.5)) +
                                         post_array[o][0],
                                         radium * (j) * math.sin(2 * 3.14159 / (j * 4) * (p + (k % 2) * 0.5)) +
                                         post_array[o][1], 0.3 * (k)])

    # brush layer
    for o in range(post_number):
        for k in range(t):
            for p in range(0, 4):
                AtomList.append([M * chain_length + Post_array_number + 4 * t * o + 4 * k + p + 1, M + 3, 1,
                                 (2 * radium + bond2_length) * math.cos(2 * 3.14159 / (4) * (p + (k % 2) * 0.5)) +
                                 post_array[o][0],
                                 (2 * radium + bond2_length) * math.sin(2 * 3.14159 / (4) * (p + (k % 2) * 0.5)) +
                                 post_array[o][1], 0.3 * (k)])

                BondList.append([(chain_length - 1) * M + 4 * t * o + 4 * k + p + 1, 2,
                                 M * chain_length + 12 * t * o + 12 * k + p * 2 + 1 + k % 2,
                                 M * chain_length + Post_array_number + 4 * t * o + 4 * k + p + 1])

    # writing headers and data into text file for simulation input
    f = open(str(M) + 'A' + str(N) + 'B' + str(chain_length - N) + '.dat', 'w')

    f.write('LAMMPS data file: Standing single chain in a box\n')
    f.write(' \n')
    f.write(str(AtomNo) + ' atoms\n')
    f.write(str(AtomTypes) + ' atom types\n')
    f.write(str(BondNo) + ' bonds\n')
    f.write(str(BondTypes) + ' bond types\n')
    f.write(' \n')  # !/usr/bin/env python3
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


N = [2]  # number of A bead
chain_length = 10
distance = 12
number = [4, 4]#number of post in both x and y dimensions
z = number[0] * distance / 2
box = [0.0, number[0] * distance , 0.0, number[1] * distance * 1.732 / 2, 0.0, z]
xsize = box[1] - box[0]
ysize = box[3] - box[2]
density = 5.0

# the position x,y of the posts
post_list = []
for j in range(number[1]):
    for i in range(number[0]):
        if (j % 2 == 0):
            post_list.append([i * distance, j * distance * 1.732 / 2])
        else:
            post_list.append([(i + 0.5) * distance, j * distance * 1.732 / 2])
post_array = np.array(post_list)
post_number = post_array.shape[0]
t = int(z / 0.3)
for n in N:
    SingleChain(int((xsize * ysize * z * density - 4 * t * post_number) / chain_length), n, z, box, chain_length,
                post_array)  
