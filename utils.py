#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 12:37:06 2021

@author: mattia
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


def grouping1(x, y, threshold=0.15):
    if len(x) != len(y):
        raise IndexError("Input arrays must be of the same length")

    # create array as [[x1,x2..xn],[y1,y2..yn]]
    data = np.vstack((x, y)).T

    # list of indexes for each array element
    elements = list(range(len(x)))

    space = []  # space contains all the clusters

    # keep going till all elements of the array have been explored
    while len(elements) != 0:
        # first elem removed from element list and added to the first cluster
        elem = elements.pop(0)
        space.append([elem])

        neighbours = [None]  # nighbour queque fake initialization

        while neighbours:
            ref = data[elem]  # coordinates of reference element

            # remove current reference element from neighbour list to avoid
            # comparing with itself
            # except for the first loop, for which the queue must be empty
            if neighbours == [None]:
                neighbours = []
            else:
                neighbours.pop(0)

            # chech with elements are within radius from current ref
            for neigh in elements:
                dist = np.linalg.norm(ref - data[neigh])
                if dist < threshold:
                    neighbours.append(neigh)

            # if neighbour list's empty, exit loop and move to next cluster
            if not neighbours:
                break

            neighbours = list(set(neighbours))  # remove duplicates
            neighbours.sort()  # order (is it useful?)

            elem = neighbours[0]  # reference for next loop

            elements.remove(elem)  # remove item from all elements
            space[-1].append(elem)  # add to the last cluster
    return space


def grouping2(x, y, threshold=0.15):
    if len(x) != len(y):
        raise IndexError("Input arrays must be of the same length")

    # create array as [[x1,x2..xn],[y1,y2..yn]]
    data = np.vstack((x, y)).T

    # list of indexes for each array element
    elements = list(range(len(x)))

    space = []  # space contains all the clusters

    # keep going till all elements of the array have been explored
    while elements:
        elem = elements[0]  # set first element
        space.append([])  # cluster waiting to be populated

        neighbours = [elem]  # start findind neighbours from elem

        while neighbours:
            neighbours = list(set(neighbours))  # remove duplicates
            neighbours.sort()  # order (is it useful?)

            elem = neighbours.pop(0)
            elements.remove(elem)  # so that we don't check for it again
            ref = data[elem]  # coordinates of reference element

            # chech with neighbours are within radius from current elem
            ## coulb be vectorialized someway... this is slow
            for neigh in elements:

                dist = np.linalg.norm(ref - data[neigh])
                if dist < threshold:
                    neighbours.append(neigh)

            space[-1].append(elem)  # add to the current (last) cluster

    return space


def grouping3(x, y, threshold=0.15, order="abscissa"):
    """
    Vectorialized version of grouping.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    threshold : TYPE, optional
        DESCRIPTION. The default is 0.15.

    Raises
    ------
    IndexError
        DESCRIPTION.

    Returns
    -------
    space : TYPE
        DESCRIPTION.

    """
    if len(x) != len(y):
        raise IndexError("Input arrays must be of the same length")

    # create array as [[x1,x2..xn],[y1,y2..yn]]
    data = np.vstack((x, y)).T

    # list of indexes for each array element
    elements = list(range(len(x)))

    space = []  # space contains all the clusters

    # keep going till all elements of the array have been explored
    while elements:
        elem = elements[0]  # set first element
        space.append([])  # cluster waiting to be populated

        neighbours = [elem]  # start findind neighbours from elem

        while neighbours:
            neighbours = list(set(neighbours))  # remove duplicates
            neighbours.sort()  # order (is it useful?)

            elem = neighbours.pop(0)
            elements.remove(elem)  # so that we don't check for it again

            # chech with neighbours are within radius from current elem
            ## coulb be vectorialized someway... this is slow

            elems_coord = data[elements]
            ref = np.expand_dims(data[elem], axis=0)

            # compute cartesian distance between ref and all points in neighbours
            distances = np.linalg.norm(ref - elems_coord, axis=1)
            new_neigh = np.where(distances < threshold)[0]  # where returns tuple
            neighbours.extend([elements[i] for i in new_neigh.tolist()])

            space[-1].append(elem)  # add to the current (last) cluster

    if order == "size-decremental":
        space.sort(key=len, reverse=True)
    return space


def grouping(x, y, x_thr=0.15, y_thr=0.15, order="abscissa", edges=(0, 1)):
    """
    Vectorialized version of grouping.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    threshold : TYPE, optional
        DESCRIPTION. The default is 0.15.

    Raises
    ------
    IndexError
        DESCRIPTION.

    Returns
    -------
    space : TYPE
        DESCRIPTION.

    """

    if len(x) != len(y):
        raise IndexError("Input arrays must be of the same length")

    # create array as [[x1,x2..xn],[y1,y2..yn]]
    data = np.vstack((x, y)).T

    # list of indexes for each array element
    elements = list(range(len(x)))

    space = []  # space contains all the clusters

    # keep going till all elements of the array have been explored
    while elements:
        elem = elements[0]  # set first element
        space.append([])  # cluster waiting to be populated

        neighbours = [elem]  # start findind neighbours from elem

        while neighbours:
            neighbours = list(set(neighbours))  # remove duplicates
            neighbours.sort()  # order (is it useful?)

            elem = neighbours.pop(0)
            elements.remove(elem)  # so that we don't check for it again

            # chech which neighbours are within radius from current elem
            dist = np.abs(data[elements] - data[elem])
            mask = np.logical_and(dist[:, 0] < x_thr, dist[:, 1] < y_thr)

            new_neigh = np.where(mask)[0]
            neighbours.extend([elements[i] for i in new_neigh.tolist()])

            space[-1].append(elem)  # add to the current (last) cluster

    # function for weighted ordering based on number of points and position
    # in reverse mode, the less points a group has the best!
    def weighted(cluster):
        n_elements = len(x)
        cluster_pos = np.average(data[cluster, 0])
        cluster_weight = len(np.array(cluster).flatten()) / n_elements
        if order == "weighted-pos reverse":
            cluster_weight = 1 - cluster_weight
        norm_pos = normalize_data(cluster_pos, edges[0], edges[1])
        return norm_pos * cluster_weight

    if order == "size-decremental":
        space.sort(key=len, reverse=True)
    if order == "weighted-pos":
        space.sort(key=weighted, reverse=True)
    if order == "weighted-pos reverse":
        space.sort(key=weighted)

    return space


def derivative(a, b=None):
    if b is None:
        return np.diff(a)

    doubled = np.repeat(a, 2)[1:-1]
    new_xrange = np.mean(doubled.reshape(-1, 2), axis=1)
    return new_xrange, np.diff(b) / np.diff(a)


def normalize_data(data, minval=None, maxval=None):
    if minval is None:
        minval = np.min(data)
        maxval = np.max(data)

    return (data - minval) / (maxval - minval)


def find_maximas(arr, ave=1, gauss=None, prominence=0.0001):
    # Find maximas of an array
    # The array can be processed with block average and gaussian smoothing
    # Maximas are found based on the prominence parameter
    n_points = len(arr)
    rounded = n_points - (n_points % ave)
    arr = arr[:rounded]

    arr = np.mean(arr.reshape(-1, ave), axis=1)

    if gauss is not None:
        arr = gaussian_filter1d(arr, gauss)

    peaks, _ = find_peaks(arr, prominence=prominence)

    # return peak position in original array
    return peaks * ave


# moving average
# ox_volt = np.convolve(ox_volt, np.ones(moving), 'valid')/moving
# ox_curr = np.convolve(ox_curr, np.ones(moving), 'valid')/moving
