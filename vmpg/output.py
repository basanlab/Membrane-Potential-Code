#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Living Matter Group - ICS-2 Juelich Research Center <j.elgeti@fz-juelich.de>
#
# Distributed under terms of the MIT license.

import datetime
import numpy as np
import meshio
import scipy.spatial as sps
from vmpg.voronoi import calculate_geometry
from vmpg.system import System
from constants import *
"""
def write_mesh(filename, points, cells, cell_data):

    lc = len(cells)

    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 4.2\nwritten by meshio v2.3.10\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS {} double\n".format(len(points)))
        f.write(' '.join(map(str, points.flat)))
        f.write('\n')
        f.write("CELLS {} {}\n".format(lc, 4*lc))
        for c0, c1, c2 in cells:
            f.write("3 {} {} {}\n".format(c0, c1, c2))
        f.write("CELL_TYPES {}\n".format(lc))
        for _ in range(lc):
            f.write("5\n")
        f.write("CELL_DATA {}\n".format(lc))
        f.write("FIELD FieldData {}\n".format(len(cell_data)))
        for name, data in cell_data.items():
            f.write("{} 1 {} double\n".format(name, lc))
            f.write(' '.join(map(str, data.flat)))
            f.write('\n')
"""

def write_dat(filename, system):
    np.savetxt(filename, system.positions)

def read_initial_state(filename):
    m = meshio.read(filename)
    cell_data = m.cell_data['triangle']

    # The files are normall in the style name.step.vtk.
    # If that fails we will understand we are in the
    # initial step


    # Create a new system with some dummy data inside.
    # We are going to overwrite it all anyway
    N = len(np.unique(m.cell_data['triangle']['id']))
    system = System(N, [1, 1])

    # Get center of mass of the cells. We happen to know that they are the last N points
    # in the mesh.points. Also, we need to remove the 3rd dimension

    p = m.points[-N:, :2]
    system.positions = p

    # We probably can optimize this loop. Each value is going to be overwritten
    # as many times as there are triangles in each cell, which is wasteful.
    for c in range(len(cell_data['id'])):
        #Find which cell we are dealing with
        id = int(cell_data['id'][c])

        system.ion_Cl[id] = cell_data['ion_Cl'][c]
        system.ion_K[id] = cell_data['ion_K'][c]
        system.ion_Na[id] = cell_data['ion_Na'][c]

        system.chanel_Cl[id] = cell_data['chanel_Cl'][c]
        system.chanel_K[id] = cell_data['chanel_K'][c]
        system.chanel_Na[id] = cell_data['chanel_Na'][c]

        system.mass[id] = cell_data['mass'][c]
        system.height[id] = cell_data['height'][c]
        system.pressure[id] = cell_data['pressure'][c]
        system.potential[id] = cell_data['potential'][c]

        system.masscharge[id] =cell_data['masscharge'][c]
        system.tension[id] = cell_data['tension'][c]
        system.vcrit[id] = cell_data['vcrit'][c]
        system.vexpansion[id] = cell_data['vexpansion'][c]
        system.Acrit[id] =  cell_data['Acrit'][c]

    return system



def write_vtu(filename_pvd, filename, system, simulation_box, step, time, debug_mode = False):

    add_to_pvd(filename_pvd, filename, time)

    areas, perimeters, _, _, vor = calculate_geometry(system, simulation_box)
    # areas, partial_perimeters, perimeters, neighbours, vor, cells

    if debug_mode:
        file_format = "vtu-ascii"
    else:
        file_format = "vtu-binary"

    #first write all point coordinates. First vertices, than cell center-of-mass
    pts = np.zeros((len(vor.vertices) + system.num_cells, 3))
    pts[:len(vor.vertices),:2] = vor.vertices
    pts[len(vor.vertices):,:2] = system.positions

    # pts = np.vstack([vor.vertices, vor.points])
    cell_id = []
    tris = []

    area_all = []

    # These are the ion concentrations inside the cell
    ion_Cl = []
    ion_K = []
    ion_Na = []

    # prefactors in ion concentration. dry mass fraction ratio of channels vs atpase time alpha
    chanel_Cl = []
    chanel_K = []
    chanel_Na = []

    mass = []
    height = []
    pressure = []
    potential = []

    masscharge = []
    tension = []
    vcrit = []
    vexpansion = []
    Acrit = []

    # write connections/triangles ...
    # TODO: optimize this loop
    for input_idx in range(len(system.positions)):
        cell = vor.regions[vor.point_region[input_idx]]#which vorregion belongs to which cell
        if -1 in cell:
            continue
        else:
            nl = len(cell)
            for idx in range(nl):
                next_idx = (idx + 1) % nl
                #connect two vertices to the central point
                tris.append([cell[idx], cell[next_idx], input_idx + len(vor.vertices)])
                cell_id.append(input_idx)
                area_all.append(areas[input_idx])

                ion_Cl.append(system.ion_Cl[input_idx])
                ion_K.append(system.ion_K[input_idx])
                ion_Na.append(system.ion_Na[input_idx])

                chanel_Cl.append(system.chanel_Cl[input_idx])
                chanel_K.append(system.chanel_K[input_idx])
                chanel_Na.append(system.chanel_Na[input_idx])

                mass.append(system.mass[input_idx])
                height.append(system.height[input_idx])
                pressure.append(system.pressure[input_idx])
                potential.append(system.potential[input_idx])

                masscharge.append(system.masscharge[input_idx])
                tension.append(system.tension[input_idx])
                vcrit.append(system.vcrit[input_idx])
                vexpansion.append(system.vexpansion[input_idx])
                Acrit.append(system.Acrit[input_idx])

    cell_id = np.array(cell_id, dtype=float)
    area_all = np.array(area_all, dtype=float)

    ion_Cl = np.array(ion_Cl, dtype=float)
    ion_K = np.array(ion_K, dtype=float)
    ion_Na = np.array(ion_Na, dtype=float)

    chanel_Cl = np.array(chanel_Cl, dtype=float)
    chanel_K = np.array(chanel_K, dtype=float)
    chanel_Na = np.array(chanel_Na, dtype=float)

    mass = np.array(mass, dtype=float)
    height = np.array(height, dtype=float)
    pressure = np.array(pressure, dtype=float)
    potential = np.array(potential, dtype=float)

    masscharge = np.array(masscharge, dtype=float)
    tension = np.array(tension, dtype=float)
    vcrit = np.array(vcrit, dtype=float)
    vexpansion = np.array(vexpansion, dtype=float)
    Acrit = np.array(Acrit, dtype=float)

    meshio.write_points_cells(filename,
                              pts,
                              cells={"triangle": np.array(tris)},
                              cell_data={"triangle": {
                                             "id": cell_id,
                                             "area": area_all,
                                             "ion_Cl": ion_Cl,
                                             "ion_K": ion_K,
                                             "ion_Na": ion_Na,
                                             "chanel_Cl": chanel_Cl,
                                             "chanel_K": chanel_K,
                                             "chanel_Na": chanel_Na,
                                             "mass": mass,
                                             "height": height,
                                             "pressure": pressure,
                                             "potential": potential,
                                             "masscharge": masscharge,
                                             "tension": tension,
                                             "vcrit": vcrit,
                                             "vexpansion": vexpansion,
                                             "Acrit": Acrit
                                         }
                                    },
                              file_format=file_format)
    #write a configfile at the same time:
    confname=filename[:-3]+'config'
    with open(confname, 'w') as file:
        file.write(f'MEDIUM_ION_Cl = {MEDIUM_ION_Cl}  #in miliMolar   \n')
        file.write(f'MEDIUM_ION_K = {MEDIUM_ION_K}   \n')
        file.write(f'MEDIUM_ION_Na = {MEDIUM_ION_Na}  \n')
        file.write(f'FARADAY = {FARADAY} # 96485/1000./1000.0 #we are using mili molar, and also mV  \n')
        file.write(f'RG = {RG} #8.314/1000 #we are using mili molar  \n')
        file.write(f'TEMP = {TEMP}  \n')
        file.write(f'DT = {DT}  \n')
        file.write(f'GAMMA_BACKGROUND = {GAMMA_BACKGROUND} #sets the pressure relaxation time.  \n')
        file.write(f'CELL_OSMOLITES = {CELL_OSMOLITES} #addtional osmolites in the cell, compared to medium. Adds to pressure  \n')
        file.write(f'SATURATION = {SATURATION} #Satturation constant of Na-K-Atpase \n')
        file.write(f'simulation_box = {simulation_box}  \n')
        file.write(f'dump_every = {dump_every}  \n')
        file.write(f'num_steps = {num_steps}  \n')
        file.write(f'output_name = "{output_name}"  \n')
        file.write(f'inputfile = "{filename}" #will restart from that file!  \n')
        file.write(f'curstep = {step}  \n')

def add_to_pvd(filename_pvd, filename, time_step):
    now = datetime.datetime.now()

    add_to_pvd.file_list.append((time_step, filename))

    with open(filename_pvd, "w") as f:
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<!--\n#This file was generated by voronoi-model on " + now.isoformat() + "\n-->\n")
        f.write("<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n   <Collection>\n")

        for time, file in add_to_pvd.file_list:
            f.write(f"    <DataSet timestep=\"{time}\" group=\"\" part=\"0\" file=\"{file}\"/>\n")

        f.write("  </Collection>\n")
        f.write( "</VTKFile>\n")

add_to_pvd.file_list = []


def safetxt(filename,system,area,perimiters):
    combo=np.vstack((system.positions[:,0],system.positions[:,1],area,perimiters,system.potential,
                    system.mass,system.chanel_Na,system.pressure,
                    system.ion_Cl,system.ion_K,system.ion_Na))
    np.savetxt(filename, combo.T)
#datfile, columns: 1=x,2=y,3=Area,4=perimiter, 5=U,6=M,7=chanelNA,8=P,9=Cl,10=K,11=Na
