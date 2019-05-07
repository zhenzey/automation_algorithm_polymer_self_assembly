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

    ptol = 0.72
    box = [-1.5, 31.5, 0.0, 30, -2.0, 8.0]
    number = 21
    difference_array = np.array([[0.0, 1.0, 1.0, 1.0, 1.0], [1.0, 0.0, 2.0, 2.0, 2.0], [1.0, 2.0, 0.0, 2.0, 2.0], [1.0, 2.0, 2.0, 0.0, 2.0], [1.0, 2.0, 2.0, 2.0, 0.0]])
    samples = 100

    # generate the period
    box_size = []
    period = []
    for i in range(3):
        box_size.append(box[2 * i + 1] - box[2 * i])
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))

    # get the desired pattern(only one pattern)
    node0 = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100_10__25", multiple_frames=True)
    data0 = node0.compute(300)
    x0 = np.array(data0.particles['Position'])
    type0 = np.array(data0.particles['Particle Type'])
    typef0 = staircase_function(x0, type0, difference_array, period, ptol, number, box)

    difference = []
    for i in range(51):
        name = "dump.A1B3_thin_film_with_wall_XN_100__" + str(i)
        node = ovito.io.import_file(name, multiple_frames=True)
        diff = difference_function(node, typef0, period, difference_array, ptol, samples, number, box)
        difference.append(diff)
        print(diff)
    print(difference)
