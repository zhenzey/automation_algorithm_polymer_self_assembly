import numpy as np
import staircase_function
import difference_function_a_d
from numpy.linalg import norm
import time
import math
import os
import ovito
from ovito.io import import_file
from staircase_function import staircase_function
from difference_function_a_d import difference_function


def lammps(typef0, period, parameter, difference_array, ptol, samples, number, box):
    """
    #run the simulation in lammmps
    :param typef0: staircase function of the desired pattern
    :type typef0: numpy.ndarray
    :param period: the list of the periodicity in 3 dimension
    :type period: list
    :param parameter: the parameters in the system
    :type parameter: numpy.ndarray(shape(1))
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :param ptol: position tolerance
    :type ptol: float
    :param sample_fre: sample frequency
    :type sample_fre：int
    :para difference0: difference between typef0 and typef0
    :type difference0: float
    :return:
    """
    #the code should be modified if the parameter is not one dimension
    var = parameter[0,0]
    os.environ['var'] = str(var)
    os.system('sed -in-place -e "s/1 2 51.25/1 2 $var/" in.asymmetric')#modify according to the parameters
    os.system("lmp_serial -in in.asymmetric")
    node = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100", multiple_frames=True)
    difference = difference_function(node, typef0, period, difference_array, ptol, samples, number, box)
    os.system('sed -in-place -e "s/1 2 $var/1 2 51.25/" in.asymmetric')#reset the parameters
    print(difference)
    return difference