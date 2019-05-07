import numpy as np
import pyswarm
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
from pyswarm import pso

def difference(x):
    """
    #calculate the object function of different particles
    :param x: array of the parameters
    :type x numpy adarray
    :return:
    """
    difference = []
    for i in range(x.shape[0]):
        difference.append(simulation(x[i]))
    return np.array(difference)



if __name__ == "__main__":

    # dpositions, dtype, initial parameter, minimiser, difference tolerance, position tolerance, update tolerance, update rate
    # list of xlow, xhi, ylo, yhi, zlo, zhi, number of setting points in every dimension, difference_array
    parameter0 = np.array([15])
    dtol = 5.0
    ptol = 0.6288#ptol should be smaller than half of the period in every dimension
    ntol = 0.1
    newr = 0.5
    box = [-1.5, 14.6, 0, 13.1, -2.0, 7.0]
    number = 11
    difference_array = np.array([[0.0, 1.0, 1.0, 1.0, 1.0], [1.0, 0.0, 2.0, 2.0, 2.0], [1.0, 2.0, 0.0, 2.0, 2.0], [1.0, 2.0, 2.0, 0.0, 2.0], [1.0, 2.0, 2.0, 2.0, 0.0]])
    sample_fre = 10
    samples = 50

    #generate the box
    box_size = []
    period = []
    for i in range(3):
        box_size.append(box[2 * i + 1] - box[2 * i])
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100_0", multiple_frames=True)
    data0 = node0.compute(200)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period, ptol, number, box)

    def simulation(parameter):
        return lammps(typef0, period, parameter, difference_array, ptol, samples, number, box)

    lb = [15]
    ub = [87.5]
    ieqcons = []
    f_ieqcons = None
    args = ()
    kwargs = {}
    swarmsize = 5
    omega = 0.5
    phip = 0.5
    phig = 0.5
    maxiter = 100
    minstep = 1.0
    minfunc = 2.0
    debug = False
    print(pso(difference, lb, ub, ieqcons, f_ieqcons, args, kwargs, swarmsize, omega, phip, phig, maxiter, minstep,
    minfunc, debug))