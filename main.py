#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Living Matter Group - ICS-2 Juelich Research Center <j.elgeti@fz-juelich.de>
#
# Distributed under terms of the MIT license.

# import matplotlib.pyplot as plt
# import numpy as np
# import scipy as sp

from constants import *
from vmpg import voronoi
from vmpg import solve_ions
from vmpg import pressure
from vmpg import euler
from vmpg import output
from vmpg import boundaries
from vmpg.system import System
from vmpg import checks
from vmpg import divisions
from vmpg import gradiente
import numpy as np
import os
import shutil
from numpy.random import default_rng

def main():
    tote=0
    rng = default_rng()
    #Hard wired changes of parameters, be carefull!
    #This actually should go into initialization, but usefull for experimenting
    #Remove when continuing an existing run
    system = output.read_initial_state(inputfile)
    system.vexpansion[:] *=0.03 #only while adapting to new units!

    
    system.tension[:] =100*np.ones(system.num_cells)
    #adding some heterogeneity
    system.mass[:] = 180*np.ones(system.num_cells)+ 20*rng.standard_normal(system.num_cells)
    system.Acrit[:] = 200*np.ones(system.num_cells)+rng.normal(size=system.num_cells,loc=200,scale=50)  #Adder model, sigma=mean/2


    #Channel parameters for
    #WT:
    system.chanel_Na[:] = 0.000275275940181926462389100374*np.ones(system.num_cells)
    system.chanel_K[:] = 0.002062484127575191575343360739*np.ones(system.num_cells)
    #Gramacidine:
    #system.chanel_Na[:] =1.2*(0.000275275940181926462389100374)*np.ones(system.num_cells)
    #system.chanel_K[:] =(0.2 *0.000275275940181926462389100374+0.002062484127575191575343360739)*np.ones(system.num_cells)
    #Valinomycine
    #system.chanel_Na[:] = 0.000275275940181926462389100374*np.ones(system.num_cells)
    #system.chanel_K[:] = 2*0.002062484127575191575343360739*np.ones(system.num_cells)
    #remember to check solveions.py


    
    system.vcrit[:] = -61.1916849*np.ones(system.num_cells)
    system.vcrit[:] = -70*np.ones(system.num_cells)

    print(curstep)
    system.masscharge = 1200*np.ones(system.num_cells)
    for step in range(curstep,num_steps):
        if (step==0) :
           facx=200/simulation_box[0]
           facy=200/simulation_box[1]
           for idx in range(system.num_cells):
                system.positions[idx][0] *=facx#1.003
                system.positions[idx][1] *=facy#1.003
           simulation_box[0] *= facx#1.003  #1.003**231 entspricht *2
           simulation_box[1] *= facy#1.003  #1.003**231 entspricht *2
        if (step>6000 and step < 1231) :
           for idx in range(system.num_cells):
                system.positions[idx][0] /=1.003
                system.positions[idx][1] /=1.003
           simulation_box[0] /= 1.003  #1.003**231 entspricht *2
           simulation_box[1] /= 1.003  #1.003**231 entspricht *2
        if step % dump_every == 0:
            #print(step)
            outfile = "outputs/{}.{:04}.vtu".format(output_name, step)
            output.write_vtu("outputs.pvd", outfile, system, simulation_box, step, step*DT)

            areas, perimeters, _, _, _ = voronoi.calculate_geometry(system, simulation_box)
            output.safetxt("outputs/{}.{:04}.dat".format('system', step),system,areas,perimeters)
            print(step,system.num_cells,np.sum(system.mass)/np.sum(areas),np.average(system.potential),np.average(system.mass),tote)
            with open('tote.dat','a')as f:
                info=str(step)+" "+str(tote)+"\n"
                info='{} {} {} {} {}\n'.format(step,system.num_cells,np.sum(system.mass)/np.sum(areas),np.average(system.potential),tote)
                f.write(info)
                f.close()


        areas, perimeters, neighbours_list, shared_lengths, _ = voronoi.calculate_geometry(system, simulation_box)
        solve_ions.calculate_ions_steady_state(system, areas, areas*system.height)
        pressure.calculate_osmotic_pressure(system)
        forces = gradiente.calcgradient(system,simulation_box)
        euler.advance_positions(system, forces)
        boundaries.periodicbc(system,simulation_box)
        areas, perimeters, _, _, _ = voronoi.calculate_geometry(system, simulation_box)
        divisions.grow(system)
        #system = divisions.divideanddiebyarea(areas,system,simulation_box)
        #system = divisions.divideanddiebyareamass(areas,system,simulation_box)
        system,tote = divisions.divideanddiebymassPot(areas,system,simulation_box,rng,tote)
        boundaries.periodicbc(system,simulation_box)




if __name__ == "__main__":
    main()
