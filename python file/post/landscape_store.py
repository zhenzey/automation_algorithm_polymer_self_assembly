import numpy as np
from numpy.linalg import norm
import time
import math
import os
import ovito
from ovito.io import import_file

def lammps(parameter,iteration):
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
    :return:
    """
    #the code should be modified if the parameter is not one dimension
    var = parameter
    os.environ['var'] = str(var)
    os.environ['iteration'] = str(iteration)
    os.system('sed -in-place -e "s/1 2 51.25/1 2 $var/" in.asymmetric')#modify according to the parameters
    os.system('sed -in-place -e "s/XN_100/XN_100__$iteration/" in.asymmetric')
    os.system("lmp_serial -in in.asymmetric")
    os.system('sed -in-place -e "s/1 2 $var/1 2 51.25/" in.asymmetric')#reset the parameters
    os.system('sed -in-place -e "s/XN_100__$iteration/XN_100/" in.asymmetric')
    return True




if __name__ == "__main__":

    for i in range(51): # XN = 0 ~ 200
        parameter = 15 + i * 1.45
        print(i, lammps(parameter, i))