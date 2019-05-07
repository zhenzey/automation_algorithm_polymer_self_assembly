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
    :param parameter: the parameters in the system
    :type parameter: numpy.ndarray(shape(1))
    :param iteration: the iteration order of current job
    :type iteration: global int
    :param density: density of system
    :type density: float
    :param chainlength: chain length of the block polymer
    :type chainlength: int
    :return:
    """
    # input parameters
    XN = parameter[0]
    thickness = parameter[1]
    M = int(thickness * 32**2 * density / chainlength)
    os.environ['iteration'] = str(iteration)
    os.environ['XN'] = str(XN)
    os.environ['thickness'] = str(thickness)
    os.environ['M'] = str(M)
    # copy the files to modify
    os.system("cp MultiBCP_with_wall_32_XN_Z.py MultiBCP_with_wall_32_XN_Z_$iteration.py")
    os.system("cp in.asymmetric in.asymmetric_$iteration")
    os.system("cp run_lammps.sh run_lammps_$iteration.sh")
    # modify relevant files
    os.system('sed -i -e "s/k = 5/k = $thickness/" MultiBCP_with_wall_32_XN_Z_$iteration.py')
    os.system("python MultiBCP_with_wall_32_XN_Z_$iteration.py")
    os.system('sed -i -e "s/3072/$M/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/1 2 51.25/1 2 $XN/" in.asymmetric_$iteration')#modify according to the parameters
    os.system('sed -i -e "s/XN_100/XN_100_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/finalasymmetric/finalasymmetric_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/XN_Z/XN_Zr_$iteration/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/xn_z.out/xn_z_$iteration.out/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/in.asymmetric/in.asymmetric_$iteration/" run_lammps_$iteration.sh')
    # run the simulation
    os.system("sbatch run_lammps_$iteration.sh")
    return None  
