import numpy as np
import staircase_function
from numpy.linalg import norm
import time
import math
import os
import ovito
from ovito.io import import_file
from staircase_function import staircase_function


def compare_type(type_array0, type_array, difference_array):
    """
    #compare two patterns, get the difference function
    #the difference between type i and type j is characterized by the difference[i][j]
    #type 0 refers to the situation that there is no particles or atoms in the position
    :param type_array0: type array for desired pattern
    :type type_array0 : numpy.ndarray
    :param type_array: type array for the pattern after simulation
    :type type_array: numpy.ndarray
    :param difference_array: an array characterizing the difference between two type of particles
    :type difference_array: numpy.ndarray
    :return:

    """
    assert (difference_array.shape[0] == difference_array.shape[1]), "The structure of the difference array is wrong"
    assert (type_array0.shape == type_array.shape), "The type array for desired pattern is incorrect"
    for i in range(difference_array.shape[0]):
        assert (difference_array[i][i] == 0.0), "The difference between same types of particles should be zero"

    differencef = []
    for i in range(type_array.shape[0]):
        differencef.append(difference_array[np.argwhere(type_array0[i] == 1)[0,0]][np.argwhere(type_array[i] == 1)[0,0]])
    return np.array(differencef)

def integration(diff):
    """
    #integrate the difference function to get the relative difference
    :param diff: the values of the difference function
    :type diff: numpy.adarray
    :return:
    """

    stair_size = 0
    l = []
    v = []
    ex = []
    for i in range(diff.shape[0] - 1):
        if np.all(diff[i + 1] - diff[i] == 0):
            stair_size += 1
        else:
            l.append(stair_size)
            # version3
            # sum up
            #v.append(diff[i])
            #ex.append((diff[i] + diff[i + 1]) / 2)
            # version4
            # sum up the square
            v.append(diff[i]**2)
            ex.append((diff[i]**2 + diff[i + 1]**2) / 2)
            stair_size = 0

    discrepency = (np.array(l) * np.array(v)).sum() + np.array(ex).sum()

    return discrepency


def samplef(node, typef0, period, difference_array, ptol, samples, number, box):
    """
    #sample the dump file to get the average difference
    :param node: objectnode
    :type node: objectnode
    :param type0: the type array of desired pattern
    :type type0: numpy.ndarray
    :param period: 3d period
    :type period: list
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float
    :param samples: the number of frames to calculate
    :type samples: int
    :type sample_fre: int
    :param number: the number of setting points in one dimension
    :type number: int
    :param box: list of xlow, xhi, ylo, yhi, zlo, zhi
    :type box: list
    :return:
    """

    frames = node.source.num_frames
    difference_sum = 0.0
    for i in range(samples): #only calculate the last 5% steps and sample according to the sample frequency
        iteration = frames - i - 1
        data = node.compute(iteration)
        positions = data.particles['Position']
        types = data.particles['Particle Type']
        typef = staircase_function(positions, types, difference_array, period, ptol, number, box)
        difference_sum += integration(compare_type(typef0, typef, difference_array))

    return difference_sum/samples
def difference_function(node, typef0, period, difference_array, ptol, samples, number, box):
    """
    #calculate the difference
    :param node: objectnode
    :type node: objectnode
    :param type0: the type array of desired pattern
    :type type0: numpy.ndarray
    :param period: 3d period
    :type period: list
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float
    :param samples: the number of frames to calculate
    :type samples: int
    :param number: the number of setting points in one dimension
    :type number: int
    :param box: list of xlow, xhi, ylo, yhi, zlo, zhi
    :type box: list
    :return:
    """

    return samplef(node, typef0, period, difference_array, ptol, samples, number, box)
