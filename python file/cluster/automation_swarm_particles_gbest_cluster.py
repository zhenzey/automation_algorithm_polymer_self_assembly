import numpy as np
import pyswarms as ps
import difference_function_a_d
import staircase_function
import lammps
import time
import math
import os
import ovito
from numpy.linalg import norm
from ovito.io import import_file
from staircase_function import staircase_function
from difference_function_a_d import difference_function
from lammps import lammps

global iterations
iterations = -1


def difference(x):
    """
    #calculate the object function of different particles
    :param x: array of the parameters
    :type x numpy adarray
    :return:
    """
    global iterations
    for i in range(x.shape[0]):
        simulation(np.array([x[i]]))
    difference = []
    for i in range(x.shape[0]):
        name = "finalasymmetric_" + str(iterations - i) + ".dat"
        while os.path.exists(name) == False:
            time.sleep(5)
        node = ovito.io.import_file(name, multiple_frames=True)
        difference.append(difference_function(node, typef0, period, difference_array, ptol, samples, number, box))
    return np.array(difference)


if __name__ == "__main__":

    # dpositions, dtype, initial parameter, minimiser, difference tolerance, position tolerance, update tolerance, update rate
    # list of xlow, xhi, ylo, yhi, zlo, zhi, number of setting points in every dimension, difference_array
    ptol = 0.72#ptol should be smaller than half of the period in every dimension
    box = [-1.5, 33.5, 0.0, 32.0, -2.0, 8.0]
    number = 21
    difference_array = np.array([[0.0, 1.0, 1.0, 1.0, 1.0], [1.0, 0.0, 2.0, 2.0, 2.0], [1.0, 2.0, 0.0, 2.0, 2.0], [1.0, 2.0, 2.0, 0.0, 2.0], [1.0, 2.0, 2.0, 2.0, 0.0]])
    samples = 50
    swarmsize = 5
    dimension = 1

    #generate the box
    box_size = []
    period = []
    for i in range(3):
        box_size.append(box[2 * i + 1] - box[2 * i])
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100_10__25", multiple_frames=True)
    data0 = node0.compute(300)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period, ptol, number, box)

    def simulation(parameter):
        global iterations
        iterations += 1
        return lammps(parameter, iterations)

    bound = (np.array([15.0]), np.array([44.0]))
    options = {'c1': 0.5, 'c2': 0.5, 'w': 0.5, 'k': 2, 'p': 2}
    optimizer = ps.single.GlobalBestPSO(swarmsize, dimension, options, bound)
    print(optimizer.optimize(difference, 100))