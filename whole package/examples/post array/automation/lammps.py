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


def lammps(parameter, iteration, numberlist, density, chainlength):
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
    var = parameter[0]
    x = var * numberlist[0]
    y = var * numberlist[1] * 1.732 / 2
    z = x / 2
    #nodeid = 55 + iteration % 5
    M = int((x * y * z * density - 4 * numberlist[0] * numberlist[1] * int(z / 0.3)) / chainlength) 
    os.environ['var'] = str(var)
    os.environ['iteration'] = str(iteration)
    os.environ['M'] = str(M)
    #os.environ['nodeid'] = str(nodeid)
    os.system("cp in.asymmetric in.asymmetric_$iteration")
    os.system("cp run_lammps.sh run_lammps_$iteration.sh")
    os.system("cp MultiBCP_hexagonal_post.py MultiBCP_hexagonal_post_$iteration.py")
    os.system('sed -i -e "s/distance = 12/distance = $var/" MultiBCP_hexagonal_post_$iteration.py')
    os.system("python MultiBCP_hexagonal_post_$iteration.py")
    time.sleep(3)
    os.system('sed -i -e "s/3072/$M/" in.asymmetric_$iteration')#modify according to the parameters
    os.system('sed -i -e "s/XN_100/XN_100_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/finalasymmetric/finalasymmetric_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/POSTA2B8/POSTA2B8_$iteration/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/posta2b8.out/posta2b8_$iteration.out/" run_lammps_$iteration.sh')
    #os.system('sed -i -e "s/node55/node$nodeid/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/in.asymmetric/in.asymmetric_$iteration/" run_lammps_$iteration.sh')
    os.system("sbatch run_lammps_$iteration.sh")
    return None  
