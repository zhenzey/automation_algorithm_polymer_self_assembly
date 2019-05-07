import numpy as np
import difference_function_d_a
import staircase_function
import time
import math
import os
import ovito
from numpy.linalg import norm
from ovito.io import import_file
from staircase_function import staircase_function
from difference_function_d_a import difference_function


if __name__ == "__main__":


    number = 31
    difference_array = np.array([[0.0, 1.0, 1.0, 1.0], [1.0, 0.0, 2.0, 2.0], [1.0, 2.0, 0.0, 2.0], [1.0, 2.0, 2.0, 0.00]])
    samples = 80
    numberlist = [4, 4]

    # generate the period
    ptol = 0.665 #ptol should be smaller than half of the period in every dimension
    box = [0.0, 48.0, 0.0, 41.568, 0.0, 24.0]
    box_size = []
    period = []
    for i in range(3):
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A2B8_thin_film_with_wall_XN_100_35", multiple_frames=True)
    data0 = node0.compute(200)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period, ptol, number, box)

    difference = []
    for i in range(51):
        distance = 5 + 0.2 * i
        ptol = distance * numberlist[1] * 1.732 / 2 * 0.48
        box = [0.0, numberlist[0] * distance, 0.0, numberlist[1] * distance * 1.732 / 2, 0.0, numberlist[0] * distance / 2]
        for i in range(3):
            period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))  
        name = "dump.A2B8_thin_film_with_wall_XN_100_" + str(i)
        node = ovito.io.import_file(name, multiple_frames=True)
        diff = difference_function(node, typef0, period, difference_array, ptol, samples, number, box)
        difference.append(diff)
        print(diff)
    print(difference)
