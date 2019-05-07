import numpy as np
import difference_function_a_d
import staircase_function
import time
import math
import os
import ovito
from numpy.linalg import norm
from ovito.io import import_file
from staircase_function import staircase_function
from difference_function_a_d import difference_function


if __name__ == "__main__":



    number = 31
    difference_array = np.array([[0.0, 1.0, 1.0, 0.0, 0.0], [1.0, 0.0, 2.0, 0.0, 0.0], [1.0, 2.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]])
    samples = 100

    # generate the period
    box0 = [-1.5, 33.5, 0.0, 32, -2.0, 8.0】
    ptol = 0.16
    box_size = []
    period0 = []
    for i in range(3):
        period0.append((box0[2 * i + 1] - box0[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100_102", multiple_frames=True)
    data0 = node0.compute(250)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period0, ptol, number, box0)

    difference = []
    for j in np.arange(5,15.5,0.5):
        for i in range(41):
            k = int(41 * (j - 5) * 2 + i)
            period = []
            box = [-1.5, 33.5, 0.0, 32.0, -2.0, 2.0 + j]
            for l in range(3):
                period.append((box[2 * l + 1] - box[2 * l]) / (number - 1))
            ptol = (j + 4) / (number - 1) * 0.48
            name = "dump.A1B3_thin_film_with_wall_XN_100_" + str(k)
            node = ovito.io.import_file(name, multiple_frames=True)
            diff = difference_function(node, typef0, period, difference_array, ptol, samples, number, box)
            difference.append(diff)
            print(diff)
    print(difference)
