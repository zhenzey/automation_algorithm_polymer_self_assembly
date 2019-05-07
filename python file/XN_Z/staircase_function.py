import numpy as np
from numpy.linalg import norm
import time
import math
import os
import ovito
from ovito.io import import_file


def staircase_function(position, particles_type, difference_array, period, ptol, number, box):
    """
    #for those particles are in the exact position of setting point, record the type, otherwise make the type 0

    :param x: the positions of particles
    :type x: numpy.ndarray
    :param particles_type: the type of particles in each position
    :type particles_type: numpy.ndarray
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float

    """
    x = np.zeros(position.shape)
    for i in range(3):
        x[:, i] = (position[:, i] - box[2 * i]) / period[i]
    a = x[:, 0]
    b = x[:, 1]
    c = x[:, 2]
    ind = np.lexsort((c, b, a))
    x = np.array([(a[i], b[i], c[i]) for i in ind])
    x_neigh = np.where((x - np.floor(x)) < 0.5, np.floor(x), np.floor(x) + 1)  # get the closest point of one bead
    particles_type = np.array([(particles_type[i]) for i in ind])

    distance = []
    distype = (np.zeros((number**3, difference_array.shape[0]))).tolist()
    distance_min = (np.ones((number**3)) * ptol).tolist()
    for i in range(x.shape[0]):
        distance.append(np.sqrt((((x[i] - x_neigh[i]) * period) ** 2).sum())) #calculate the distance to the nearest setting point
        order = int(number**2 * x_neigh[i][0] + number * x_neigh[i][1] + x_neigh[i][2])
        if x_neigh[i][2] < number and distance[i] < ptol and distance[i] < distance_min[order]:
            distype[order] = np.zeros((difference_array.shape[0])).tolist()
            distype[order][particles_type[i]] = 1
            distance_min[order] = distance[i]
    distype = np.array(distype)
    for i in range(number**3):
        if np.all(distype[i] == 0.0):
            distype[i][0] = 1

    return distype

