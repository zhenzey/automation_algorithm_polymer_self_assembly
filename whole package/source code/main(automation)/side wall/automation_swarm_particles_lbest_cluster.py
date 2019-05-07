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
        simulation(np.array(x[i]))
        time.sleep(1)
    difference = []
    for i in np.arange(x.shape[0] - 1, -1, -1):
        fname = "finalasymmetric_" + str(iterations - i) + ".dat"
        # check whether the simulation is finished
        while os.path.exists(fname) == False:
            time.sleep(5)
        name = "dump.A2B8_thin_film_with_wall_XN_100_10_" + str(iterations - i)
        node = ovito.io.import_file(name, multiple_frames=True)
        difference.append(difference_function(node, typef0, period, difference_array, ptol, samples, number, box))
    print(difference)
    return np.array(difference)


if __name__ == "__main__":


    ptol = 0.72
    box = [-1.5, 33.5, 0.0, 32.0, -2.0, 8.0]
    number = 21
    difference_array = np.array([[0.0, 1.0, 1.0, 1.0, 1.0], [1.0, 0.0, 2.0, 2.0, 2.0], [1.0, 2.0, 0.0, 2.0, 2.0], [1.0, 2.0, 2.0, 0.0, 2.0], [1.0, 2.0, 2.0, 2.0, 0.0]])
    samples = 50
    swarmsize = 5
    dimension = 1
    ftol = 1e-2
    #generate the box
    box_size = []
    period = []
    for i in range(3):
        box_size.append(box[2 * i + 1] - box[2 * i])
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A2B8_thin_film_with_wall_XN_100_10", multiple_frames=True)
    data0 = node0.compute(300)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period, ptol, number, box)

    def simulation(parameter):
        global iterations
        iterations += 1
        return lammps(parameter, iterations)


    # please refer to pyswarms module
    bound = (np.array([15.0]), np.array([44.0]))
    options = {'c1': 0.5, 'c2': 0.5, 'w': 0.9, 'p': 2, 'k': 2}
    optimizer = ps.single.LocalBestPSO(swarmsize, dimension, options, bound, None, 1, ftol)
    print(optimizer.optimize(difference, 100))
