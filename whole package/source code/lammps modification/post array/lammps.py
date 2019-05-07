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
    # copy the files to modify
    os.system("cp in.asymmetric in.asymmetric_$iteration")
    os.system("cp run_lammps.sh run_lammps_$iteration.sh")
    os.system("cp MultiBCP_hexagonal_post.py MultiBCP_hexagonal_post_$iteration.py")
    # modify relevant files
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
    # run the simulation
    os.system("sbatch run_lammps_$iteration.sh")
    return None  
