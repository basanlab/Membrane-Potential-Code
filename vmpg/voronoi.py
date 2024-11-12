#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Living Matter Group - ICS-2 Juelich Research Center <j.elgeti@fz-juelich.de>
#
# Distributed under terms of the MIT license.

import numpy as np
import scipy.spatial as sps
from vmpg import boundaries

def calculate_geometry(system, simulation_box):

    N = len(system.positions)
    periodic_positions = np.zeros((9*N, 2))

    #Main image
    periodic_positions[:N, :] = system.positions

    #"Right" side
    periodic_positions[1*N:2*N, :] = system.positions + np.array([simulation_box[0], -simulation_box[1]])
    periodic_positions[2*N:3*N, :] = system.positions + np.array([simulation_box[0], 0])
    periodic_positions[3*N:4*N, :] = system.positions + np.array([simulation_box[0], +simulation_box[1]])

    #"Middle"
    periodic_positions[4*N:5*N, :] = system.positions + np.array([0, -simulation_box[1]])
    periodic_positions[5*N:6*N, :] = system.positions + np.array([0, simulation_box[1]])

    #"Left" side
    periodic_positions[6*N:7*N, :] = system.positions + np.array([-simulation_box[0], -simulation_box[1]])
    periodic_positions[7*N:8*N, :] = system.positions + np.array([-simulation_box[0], 0])
    periodic_positions[8*N:9*N, :] = system.positions + np.array([-simulation_box[0], +simulation_box[1]])

    # Create the Voronoi tesselation. The qhull options are the default recommended by scipy.
    # According to the documentation, these mean:
    # * Qbb: scale the last coordinate to [0,m] for Delaunay
    # * Qc: keep coplanar points with nearest facet
    # * Qx: exact pre-merges (allows coplanar facets)
    vor = sps.Voronoi(periodic_positions, qhull_options='Qbb Qc Qx')

    # The cell ids returned by sps.Voronoi are not automatically related to the
    # ids of the points. The next loop fixes that
    cells = []
    for input_idx in range(N):
        cells.append(vor.regions[vor.point_region[input_idx]])

    areas = calculate_areas(vor.vertices, cells)
    partial_perimeters = calculate_partial_perimeters(vor.vertices, cells)
    perimeters = calculate_perimeters(partial_perimeters)
    neighbours, shared_lengths = calculate_neighbour_list(vor, N)

    return areas, perimeters, neighbours, shared_lengths, vor

def calculate_areas(vertices, cells):
    areas = np.zeros(len(cells))

    for idx, cell in enumerate(cells):
        # Some cells are connected to a point at infinity. This is signaled
        # by an index of -1 in the vertex list. If we have one of those, we
        # set the area as -1 for easier filtering later.
        if -1 in cell:
            areas[idx] = -1
        else:
            loc_v = vertices[cell, :]
            x = loc_v[:,0]
            y = loc_v[:,1]

            # Use the shoelace formula to calculate the area of the
            # cell, ??? added factor two???
            area = np.dot(x, np.roll(y,1)) - np.dot(y, np.roll(x,1))
            areas[idx] = np.abs(area)/2

    return areas

def calculate_partial_perimeters(vertices, cells):
    # Some cells are connected to a point at infinity. This is signaled
    # by an index of -1 in the vertex list. If we have one of those, we
    # set the perimeter as -1 for easier filtering later.
    partial_perimeters = []
    for cell in cells:
        if -1 in cell:
            partial_perimeters.append(np.ones(len(cell))*-1)
        else:
            loc_v = vertices[cell, :]
            d = loc_v - np.roll(loc_v, 1, axis=0)
            partial_perimeters.append(np.linalg.norm(d, axis=1))

    return partial_perimeters

def calculate_perimeters(partial_perimeters):
    return np.array(
        [np.sum(pp) for pp in partial_perimeters]
    )

def calculate_neighbour_list(vor, N):

    """
    Calculate the neighbor list of all cells. This will be a list of lists.
    Each edge of each voronoi cell is called a "ridge".

    The `ridge` parameter is a list where every element is a list with two
    indices. These indices are related to the index of the center points of the
    Voronoi cell. For instance, an element [3, 8] would mean that the cell
    formed by points number 3 and number 8 share a border.

    The `ridge_vertices` parameter contain the indices of the Voronoi
    vertices that are connected by the edge.

    TODO: Fix double connections
    """

    neighbours = []
    shared_lengths = []
    for i in range(N):
        neighbours.append([])
        shared_lengths.append([])

    for ridge_idx, ridge_points in enumerate(vor.ridge_points):
        # If the index of the neighbor is larger
        # than the number of particles, we are
        # connected to the mirror image. Fix
        # the index to connect to the cell in
        # the original system image.
        if(ridge_points[0]>=N and ridge_points[1]>=N):
            continue
        idx0 = ridge_points[0]%N
        idx1 = ridge_points[1]%N

        if idx0 == idx1:
            continue #should not happen?
        # Now, calculate the length of the edge
        # connecting cells idx0 and idx1. First
        # we get their indices
        v0_idx = vor.ridge_vertices[ridge_idx][0]
        v1_idx = vor.ridge_vertices[ridge_idx][1]

        # And now the actual position of the vertices
        # creating the edge
        v0 = vor.vertices[v0_idx]
        v1 = vor.vertices[v1_idx]
        # And now we calculate the length of the connecting
        # edge
        d = np.linalg.norm(v0-v1)

        # TODO: Check if we need to check this at all, and if
        # so, if searching just one of the two indices is
        # enough
        if idx1 not in neighbours[idx0] and idx0 not in neighbours[idx1]:
            neighbours[idx0].append(idx1)
            neighbours[idx1].append(idx0)
            shared_lengths[idx0].append(d)
            shared_lengths[idx1].append(d)

    #sanitiy check - should be no self neighbors
    for idx, n in enumerate(neighbours):
        if idx in n:
            raise ValueError(f"self neighbor {idx}")

    return neighbours, shared_lengths

def calculate_voronoi_energy(pos, system, simulation_box):

    N = len(pos)
    periodic_positions = np.zeros((9*N, 2))

    #Main image
    periodic_positions[:N, :] = pos

    #"Right" side
    periodic_positions[1*N:2*N, :] = pos + np.array([simulation_box[0], -simulation_box[1]])
    periodic_positions[2*N:3*N, :] = pos + np.array([simulation_box[0], 0])
    periodic_positions[3*N:4*N, :] = pos + np.array([simulation_box[0], +simulation_box[1]])

    #"Middle"
    periodic_positions[4*N:5*N, :] = pos + np.array([0, -simulation_box[1]])
    periodic_positions[5*N:6*N, :] = pos + np.array([0, simulation_box[1]])

    #"Left" side
    periodic_positions[6*N:7*N, :] = pos + np.array([-simulation_box[0], -simulation_box[1]])
    periodic_positions[7*N:8*N, :] = pos + np.array([-simulation_box[0], 0])
    periodic_positions[8*N:9*N, :] = pos + np.array([-simulation_box[0], +simulation_box[1]])

    # Create the Voronoi tesselation. The qhull options are the default recommended by scipy.
    # According to the documentation, these mean:
    # * Qbb: scale the last coordinate to [0,m] for Delaunay
    # * Qc: keep coplanar points with nearest facet
    # * Qx: exact pre-merges (allows coplanar facets)
    vor = sps.Voronoi(periodic_positions, qhull_options='Qbb Qc Qx')

    # The cell ids returned by sps.Voronoi are not automatically related to the
    # ids of the points. The next loop fixes that
    cells = []
    for input_idx in range(N):
        cells.append(vor.regions[vor.point_region[input_idx]])

    areas = calculate_areas(vor.vertices, cells)
    pzero=3.81*np.sqrt(areas)
    partial_perimeters = calculate_partial_perimeters(vor.vertices, cells)
    perimeters = calculate_perimeters(partial_perimeters)
    energy=0
    for input_idx in range(N):
        energy+= -system.pressure[input_idx]*system.height[input_idx]*areas[input_idx]\
                 +system.tension[input_idx]*(perimeters[input_idx]-pzero[input_idx])**2


    return energy
