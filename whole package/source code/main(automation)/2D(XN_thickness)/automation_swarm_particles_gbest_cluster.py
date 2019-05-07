import numpy as np
import pyswarmsmod as ps
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
    :param x: array of the parameters, standard form x = [100(XN), 6(thickness)]
    :type x numpy adarray
    :return:
    """
    global iterations
    for i in range(x.shape[0]):
        simulation(np.array(x[i]))
        time.sleep(1)
    difference = []
    ptol = 0.512
    for i in np.arange(x.shape[0] - 1, -1, -1):
        fname = "finalasymmetric_" + str(iterations - i) + ".dat"
        # check whether the simulation is finished
        while os.path.exists(fname) == False:
            time.sleep(5)
        name = "dump.A2B8_thin_film_with_wall_XN_100_" + str(iterations - i)
        node = ovito.io.import_file(name, multiple_frames=True)
        difference.append(difference_function(node, typef0, period, difference_array, ptol, samples, number, box))
    print(difference)
    return np.array(difference)


if __name__ == "__main__":


    density = 5
    chainlength = 10
    number = 31
    difference_array = np.array([[0.0, 1.0, 1.0, 0.0, 0.0], [1.0, 0.0, 2.0, 0.0, 0.0], [1.0, 2.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]])
    samples = 100
    swarmsize = 5
    dimension = 2
    ftol = 1e-2 # tolerance to jump out of the iteration

    #generate the box
    period = []
    box = [-1.5, 33.5, 0.0, 32.0, -2.0, 5.0]
    ptol = 0.512
    for i in range(3):
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A2B8_thin_film_with_wall_XN_100_reference", multiple_frames=True)
    data0 = node0.compute(250)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period, ptol, number, box)

    def simulation(parameter):
        global iterations
        iterations += 1
        return lammps(parameter, iterations, density, chainlength)

    # please refer to pyswarms module
    bound = (np.array([15.0, 5.0]), np.array([44.0, 15.0]))
    options = {'c1': 0.5, 'c2': 0.5, 'w': 0.9}
    optimizer = ps.single.GlobalBestPSO(swarmsize, dimension, options, bound, None, 1, ftol)
    print(optimizer.optimize(difference, 50))
