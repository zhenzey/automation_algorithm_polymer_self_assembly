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


def lammps(parameter, iteration):
    """
    #run the simulation in lammmps
    :param parameter: the parameters in the system
    :type parameter: numpy.ndarray(shape(1))
    :param iteration: the iteration order of current job
    :type iteration: global int
    """
    # input parameters
    var = parameter[0]
    os.environ['var'] = str(var)
    os.environ['iteration'] = str(iteration)
    # copy the files to modify
    os.system("cp in.asymmetric in.asymmetric_$iteration")
    os.system("cp run_lammps.sh run_lammps_$iteration.sh")
    # modify the relevant files
    os.system('sed -i -e "s/1 2 51.25/1 2 $var/" in.asymmetric_$iteration')#modify according to the parameters
    os.system('sed -i -e "s/XN_100_10/XN_100_10_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/finalasymmetric/finalasymmetric_$iteration/" in.asymmetric_$iteration')
    os.system('sed -i -e "s/A2B8singlelayer/A2B8singlelayer_$iteration/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/a2b8singlelayer.out/a2b8singlelayer_$iteration.out/" run_lammps_$iteration.sh')
    os.system('sed -i -e "s/in.asymmetric/in.asymmetric_$iteration/" run_lammps_$iteration.sh')
    # run the simulation
    os.system("sbatch run_lammps_$iteration.sh")
    return None  
