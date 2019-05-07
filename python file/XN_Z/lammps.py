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


def lammps(parameter, iteration, density, chainlength):
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
    XN = parameter[0]
    thickness = parameter[1]
    M = int(thickness * 32**2 * density / chainlength)
    os.environ['iteration'] = str(iteration)
    os.environ['XN'] = str(XN)
    os.environ['thickness'] = str(thickness)
    os.environ['M'] = str(M)
    os.system("cp MultiBCP_with_wall_32_XN_Z.py MultiBCP_with_wall_32_XN_Z_$iteration.py")
    os.system("cp in.asymmetric in.asymmetric_$iteration")
    os.system("cp run_lammps.sh run_lammps_$iteration.sh")
    os.system('sed -i -e "s/k = 5/k = $thickness/" MultiBCP_with_wall_32_XN_Z_$iteration.py')
    os.system("python MultiBCP_with_wall_32_XN_Z_$iteration.py")
    os.system('sed -i -e "s/3072/$M/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/1 2 51.25/1 2 $XN/" in.asymmetric_$iteration')#modify according to the parameters
    os.system('sed -i -e "s/XN_100/XN_100_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/finalasymmetric/finalasymmetric_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/XN_Z/XN_Zr_$iteration/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/xn_z.out/xn_z_$iteration.out/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/in.asymmetric/in.asymmetric_$iteration/" run_lammps_$iteration.sh')
    os.system("sbatch run_lammps_$iteration.sh")
    return None  
