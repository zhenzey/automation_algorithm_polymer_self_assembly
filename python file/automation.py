import numpy as np
from numpy.linalg import norm
import time
import math
import os
import ovito
from ovito.io import import_file



def staircase_function(position, particles_type, difference_array, period, ptol, number, box):
    """
    #for those particles are in the exact position of setting point, record the type, otherwise make the type 0

    :param x: the positions of particles
    :type x: numpy.ndarray
    :param particles_type: the type of particles in each position
    :type particles_type: numpy.ndarray
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float

    """
    x = np.zeros(position.shape)
    for i in range(3):
        x[:, i] = (position[:, i] - box[2 * i]) / period[i]
    a = x[:, 0]
    b = x[:, 1]
    c = x[:, 2]
    ind = np.lexsort((c, b, a))
    x = np.array([(a[i], b[i], c[i]) for i in ind])
    x_neigh = np.where((x - np.floor(x)) < 0.5, np.floor(x), np.floor(x) + 1)  # get the closest point of one bead
    particles_type = np.array([(particles_type[i]) for i in ind])

    distance = []
    distype = (np.zeros((number**3, difference_array.shape[0]))).tolist()
    distance_min = (np.ones((number**3)) * ptol).tolist()
    for i in range(x.shape[0]):
        distance.append(np.sqrt((((x[i] - x_neigh[i]) * period) ** 2).sum())) #calculate teh distance to the nearest setting point
        order = int(number**2 * x_neigh[i][0] + number * x_neigh[i][1] + x_neigh[i][2])
        if distance[i] < ptol and distance[i] < distance_min[order]:
            distype[order] = np.zeros((difference_array.shape[0])).tolist()
            distype[order][particles_type[i]] = 1
            distance_min[order] = distance[i]
    distype = np.array(distype)
    for i in range(number**3):
        if np.all(distype[i] == 0.0):
            distype[i][0] = 1

    return distype




   # set_x = []
    #type = []
    #interval = ptol / np.array(period)# largest interval in 3 dimensions
    #for i in range(number):#z axis direction
    #   for j in range(number): #y axis direction
    #        for k in range(number): #x axis direction
    #           # set_x.append([k * period[0] + box[0], j * period[1] + box[1], i * period[2] + box[2]])
    #           l = 0
    #            distance = []
    #            while (l < x.shape[0] and np.all((x[l, :] - np.array([k, j, i]) - interval) < 0)):
    #                distance.append(math.sqrt((((x[l, :] - np.array([k, j, i])) * period) ** 2).sum()))
    #                l += 1
    #                if l < x.shape[0] and np.all((x[l, :] - np.array([k, j, i]) - interval) > 0):
    #                    break
    #            if len(distance) != 0 and distance[np.argmin(distance, axis=0)] < ptol:
    #               type.append(particles_type[np.argmin(distance, axis=0)])
    #            else:
    #               type.append(0.0)
    #return np.array(type)




    # modify configurations to find the closest setting point to a bead
    #x = np.zeros(position.shape)
    #for i in range(3):
    #   x[:, i] = position[:, i] / period[i]
    #a = x[:,0]
    #b = x[:,1]
    #c = x[:,2]
    #ind = np.lexsort((c, b, a))
    #x = np.array([(a[i], b[i], c[i]) for i in ind])
    #particles_type = np.array([(particles_type[i]) for i in ind])

    #x_neigh = np.where((x - np.floor(x)) < 0.5, np.floor(x), np.floor(x) + 1)  # get the closest point of one bead

    # decide whether the closest point is overlapped by the bead
    #distance = np.sqrt((((x - x_neigh) * period) ** 2).sum(1))  # calculate every bead's distance to its nearest point
    #distance_chain = np.where(distance < ptol, 1, 0)

    #return distance_chain * particles_type


def compare_type(type_array0, type_array, difference_array):
    """
    #compare two patterns, get the difference function
    #the difference between type i and type j is characterized by the difference[i][j]
    #type 0 refers to the situation that there is no particles or atoms in the position
    :param type_array0: type array for desired pattern
    :type type_array0 : numpy.ndarray
    :param type_array: type array for the pattern after simulation
    :type type_array: numpy.ndarray
    :param difference_array: an array characterizing the difference between two type of particles
    :type difference_array: numpy.ndarray
    :return:

    """
    assert (difference_array.shape[0] == difference_array.shape[1]), "The structure of the difference array is wrong"
    assert (type_array0.shape == type_array.shape), "The type array for desired pattern is incorrect"
    for i in range(difference_array.shape[0]):
        assert (difference_array[i][i] == 0.0), "The difference between same types of particles should be zero"
    type_array_shift = np.zeros(type_array.shape)
    difference = 0.0
    for i in range(difference_array.shape[0]):
        difference_array_shift = []
        for j in range(difference_array.shape[0]):
            type_array_shift[:, (j + i) % 5] = type_array[:, j]
            difference_array_shift.append(difference_array[j][(j + i) % 5])
        difference += (type_array_shift * type_array0 * np.array(difference_array_shift)).sum()
    return difference


def integration(diff):
    """
    #integrate the difference function to get the relative difference
    :param diff: the values of the difference function
    :type diff: numpy.adarray
    :return:
    """

    stair_size = 0
    l = []
    v = []
    ex = []
    for i in range(diff.shape[0] - 1):
        if diff[i + 1] - diff[i] == 0:
            stair_size += 1
        else:
            l.append(stair_size)
            v.append(diff[i])
            ex.append((diff[i] + diff[i + 1]) / 2)
            stair_size = 0

    discrepency = (np.array(l) * np.array(v)).sum() + np.array(ex).sum()

    return discrepency


def samplef(node, period, difference_array, ptol, samples, number, box):
    """
    #sample the dump file to get the average difference
    :param node: objectnode
    :param type0: the type array of desired pattern
    :type type0: numpy.ndarray
    :param period: 3d period
    :type period: list
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarraynode, period, difference_array, ptol, samples, number, box
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float
    :param sample_fre: sample frequency
    :type sample_fre: int
    :return:
    """

    frames = node.source.num_frames
    typef = np.zeros((number ** 3, difference_array.shape[0]))
    for i in range(samples): #only calculate the last 5% steps and sample according to the sample frequency
        iteration = frames - i - 1
        data = node.compute(iteration)
        print(iteration)
        positions = data.particles['Position']
        types = data.particles['Particle Type']
        typef += staircase_function(positions, types, difference_array, period, ptol, number, box)


    return typef/samples

def lammps(typef0, period, parameter, difference_array, ptol, samples, number, box, difference0):
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
    var = parameter[0]
    os.environ['var'] = str(var)
    os.system('sed -in-place -e "s/1 2 51.25/1 2 $var/" in.asymmetric')#modify according to the parameters
    os.system("lmp_serial -in in.asymmetric")
    node = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100", multiple_frames=True)
    typef = samplef(node, period, difference_array, ptol, samples, number, box)
    print(np.all(typef == typef0))
    difference = compare_type(typef0, typef, difference_array) - difference0
    os.system('sed -in-place -e "s/1 2 $var/1 2 51.25/" in.asymmetric')#reset the parameters
    print(difference)
    return difference


def gradient(parameter0, typef0, difference0, period, difference_array, ptol, samples, delta, number):
    """
    #calculate gradient
    :param parameter0: the parameters(variables) for the point to calculate
    :type parameter0: numpy.ndarray
    :param typef0: staircase function of the desired pattern
    :type typef0: numpy.ndarray
    :param period: the list of the periodicity in 3 dimension
    :type period: list
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :param ptol: position tolerance
    :type ptol: float
    :param sample_fre: sample frequency
    :type sample_fre：int
    :param delta: the delta parameters in the differentiation
    :type delta: float
    :return:
    """
    gradient = []
    parameter_delta = []
    diff0 = lammps(typef0, period, parameter0, difference_array, ptol, samples, number, box, difference0)
    print("gradient 0")
    print(parameter0)

    for i in range(len(parameter0.tolist())): #if the shape of parameter0 is not(,1), it should be modified
        parameter_delta.append(parameter0[i] + delta)
        difference = lammps(typef0, period, np.array([parameter_delta[i]]), difference_array, ptol, samples, number,
                            box, difference0)
        print("gradient")
        gradient.append((difference - diff0) / delta)

    return np.array(gradient)


def line_search(line_f, diff0, alpha0, m, c=0.5, t=0.9):
    """

    :param line_f: in our case, it is the difference function
    :type line_f: callable function
    :param alpha: parameter for the line search
    :type alpha: float
    :param m: gradient of the difference function
    :type m:
    :param c: parameters of backtracking line search
    :type c: float
    :param t: parameters of backtracking line search
    :type t: float
    :return:
    """
    while line_f(alpha0) > diff0 - alpha0 * c * m:
        alpha0 *= t
    return alpha0

def steepest_descent(parameter0, typef0, difference0, df, period, difference_array, samples, ptol, delta, number, iterations, tolerance, c, t):
    """
    # minimize and find the local minimum
    :param parameter0: the parameters(variables) for the point to calculate
    :type parameter0: numpy.ndarray
    :param df: gradient function
    :type df: callable function
    :param period: the list of the periodicity in 3 dimension
    :type period: list
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :param ptol: position tolerance
    :type ptol: float
    :param sample_fre: sample frequency
    :type sample_fre：int
    :param delta: the delta parameters in the differentiation
    :type delta: float
    :param iterations: iterations for steepest descent
    :type iterations: int
    :param tolerance: tolerance for distance
    :type tolerance: float
    :param c: parameters of backtracking line search
    :type c: float
    :param t: parameters of backtracking line search
    :type t: float
    :return:
    """

    distance = 1e10
    alpha = 20.0
    parameter1 = parameter0.copy()
    iteration = 0


    while distance > tolerance:
        if iteration > iterations:
            return parameter1, False
        iteration += 1
        parameter0[:] = parameter1
        gradient = df(parameter0, typef0, difference0, period, difference_array, ptol, samples, delta, number)
        m = np.linalg.norm(gradient)
        p = -gradient / m
        print(gradient,p)

        def line(alpha):
            parameter = parameter0.copy()
            parameter[:] = parameter0 + alpha * p
            difference = lammps(typef0, period, parameter, difference_array, ptol, samples, number, box, difference0)
            print("line")
            return difference

        diff0 = lammps(typef0, period, parameter0, difference_array, ptol, samples, number, box, difference0)
        alpha = line_search(line, diff0, alpha, m, c, t)
        parameter1[:] = parameter0 + alpha * p
        distance = np.linalg.norm(parameter1 - parameter0)

    return parameter1, True

def minimiserf(parameter, typef0, difference0, gradient, period, difference_array,samples, ptol, number):
    return steepest_descent(parameter, typef0, difference0, gradient, period, difference_array, samples, ptol, 2.0, number, 5, 1.0, 0.1,
                            0.9)[0]


def automation(node0, parameter0, minimiserf, dtol, ptol, ntol, newr, box, number, difference_array, sample_fre, samples):
    """
    main loop
    :param x0: the positions of particles in desired pattern
    :type x0: numpy.ndarray
    :param type0: the type of particles in desired pattern
    :type type0: numpy.ndarray
    :param parameter0: the initial parameters array in the system
    :type parameter0: numpy.array
    :param dtol: tolerance for the value of difference function
    :type dtol: float
    :param ptol: postion tolerance: determining whether two positions are the same
    :type ptol: float
    :param ntol: new tolerance
    :type ntol: float
    :param newr: update rate
    :type newr: float
    :param box: list of xlo, xhi, ylo, yhi, zlo, zhi
    :type box: list
    :param number: the number of partition in every dimension
    :type number: int
    :param difference_array: the array for pairwise difference
    :type difference_array: numpy.ndarray
    :return:
    """
    # generate setting points according to the box size
    assert (number > 1), "the number of points is improper"
    step_size = 5.0
    temperature = 1.0
    new = 0
    accept = 0
    new_flag = False
    iteration = 0

    #generate the setting point
    box_size = []
    period = []
    for i in range(3):
        box_size.append(box[2 * i + 1] - box[2 * i])
        period.append((box[2 * i + 1] - box[2 * i]) / (number - 1))



    # run the first simulation
    typef0 = samplef(node0, period, difference_array, ptol, samples, number, box)#the type array for the desired pattern
    difference0 = compare_type(typef0, typef0, difference_array)
    print(difference0)
    difference = lammps(typef0, period, parameter0, difference_array, ptol, samples, number, box, difference0)
    print("first")
    difference_g_min = difference
    parameter_g_min = parameter0


    # main loop
    # smart monte carlo simulation
    while (difference > dtol):
        parameter = parameter0 + step_size * (0.5 - np.random.rand(*parameter0.shape))
        parameter_min = minimiserf(parameter, typef0, difference0, gradient, period, difference_array, samples, ptol, number)[0]
        difference_new = lammps(typef0, period, parameter_min, difference_array, ptol, samples, number, box, difference0)
        print("new")
        if abs(difference - difference_new) > ntol:
            new += 1
            new_flag = True
        else:
            new_flag = False
        if difference_new < difference or np.exp(-(difference_new - difference)/temperature) > np.random.rand(): #metropolis
            parameter0[:] = parameter_min
            difference = difference_new
            if new_flag:
                accept += 1
            if difference < difference_g_min:
                difference_g_min = difference
                parameter_g_min[:] = parameter0[:]
            # update temperature and step_size according to the accept rate and new rate
            if iteration % sample_fre == 1:
                if new / iteration > newr:
                    step_size /= 1.01
                else:
                    step_size *= 1.01
                #if accept / new > acceptr:
                 #   temperature /= 1.01
                #else:
                 #   temperature *= 1.01
            iteration += 1

    return iteration, difference_g_min, parameter_g_min


if __name__ == "__main__":
    # dpositions, dtype, initial parameter, minimiser, difference tolerance, position tolerance, update tolerance, update rate
    # list of xlow, xhi, ylo, yhi, zlo, zhi, number of setting points in every dimension, difference_array
    node0 = ovito.io.import_file("dump.A1B3_thin_film_with_wall_XN_100_0", multiple_frames=True)
    #node.compute(40)
    #x0 = node.output.particle_properties.position.array
    #type0 = node.output.particle_properties.particle_type.array
    parameter0 = np.array([15])
    dtol = 5.0
    ptol = 0.6288#ptol should be smaller than half of the period in every dimension
    ntol = 0.1
    newr = 0.5
    box = [-1.5, 14.6, 0, 13.1, -2.0, 7.0]
    number = 11
    difference_array = np.array([[0.0, 1.0, 1.0, 1.0, 1.0], [1.0, 0.0, 2.0, 2.0, 2.0], [1.0, 2.0, 0.0, 2.0, 2.0], [1.0, 2.0, 2.0, 0.0, 2.0], [1.0, 2.0, 2.0, 2.0, 0.0]])
    sample_fre = 10
    samples = 50
    print("node0, parameter0, minimiserf, dtol, ptol, ntol, newr, box, number, difference_array, sample_frequency\n")
    print(automation(node0, parameter0, minimiserf, dtol, ptol, ntol, newr, box, number, difference_array, sample_fre, samples))














