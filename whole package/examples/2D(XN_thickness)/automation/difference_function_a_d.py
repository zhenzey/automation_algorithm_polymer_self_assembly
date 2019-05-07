import numpy as np
import staircase_function
from numpy.linalg import norm
import time
import math
import os
import ovito
from ovito.io import import_file
from staircase_function import staircase_function

def samplef(node, period, difference_array, ptol, samples, number, box):
    """
    #sample the dump file to get the average difference
    :param node: objectnode
    :type node: objectnode
    :param period: 3d period
    :type period: list
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float
    :param samples: the numder of frames to calculate
    :type samples: int
    :param number: the number of setting points in one dimension
    :type number: int
    :param box: list of xlow, xhi, ylo, yhi, zlo, zhi
    :type box: list
    :return:
    """

    frames = node.source.num_frames
    typef = np.zeros((number ** 3, difference_array.shape[0]))
    for i in range(samples): #only calculate the last 5% steps and sample according to the sample frequency
        iteration = frames - i - 1
        data = node.compute(iteration)
        positions = data.particles['Position']
        types = data.particles['Particle Type']
        typef += staircase_function(positions, types, difference_array, period, ptol, number, box)

    return typef/samples

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
    # version1
    # multiply the histogram
    type_array_shift = np.zeros(type_array.shape)
    difference = 0.0
    for i in range(difference_array.shape[0]):
        difference_array_shift = []
        for j in range(difference_array.shape[0]):
            type_array_shift[:, (j + i) % 5] = type_array[:, j]
            difference_array_shift.append(difference_array[j][(j + i) % 5])
        difference += (type_array_shift * type_array0 * np.array(difference_array_shift)).sum()
    return difference
    # version2
    # chi-square
    #return np.where(type_array + type_array0 == 0, 0.0, (type_array - type_array0)**2/(type_array0 + type_array)).sum()

def difference_function(node, typef0, period, difference_array, ptol, samples, number, box):
    """
    #calculate the difference
    :param type0: the type array of desired pattern
    :type type0: numpy.ndarray
    :param node: objectnode
    :type node: objectnode
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
    return compare_type(typef0, samplef(node, period, difference_array, ptol, samples, number, box), difference_array)

