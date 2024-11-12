#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Living Matter Group - ICS-2 Juelich Research Center <j.elgeti@fz-juelich.de>
#
# Distributed under terms of the MIT license.

import numpy as np

class System:

    """Main class, holding numpy arrays for all the cells in the system."""

    def __init__(self, N, L):
        self.num_cells = N #number of cells
        self.positions =  np.random.rand(N,2)  #position of cells, micrometers
#        self.positions = np.array([[5.,5.],[15.,5.],[25.,5.],
#                                   [5.,15.],[19,15.],[25.1,15.],
#                                   [5.1,25.4],[15.1,25.],[25.1,25.]
#                                   ]) #np.random.rand(N,2)*L  #position of cells, micrometers




        # These are the ion concentrations inside the cell
        self.ion_Cl = 13.45*np.ones(N)
        self.ion_K = 150*np.ones(N)
        self.ion_Na = 10*np.ones(N)

        # prefactors in ion concentration. dry mass fraction ratio of channels vs atpase time alpha
        self.chanel_Cl = 0.006*np.ones(N)
        self.chanel_K = 0.0005499052513181528*np.ones(N)
        self.chanel_Na = 8.569452513649223e-05*np.ones(N)

        self.mass = 0.2*1000*np.ones(N)  #currently total dry mass, switched to % *Volume, will adapt
        self.height = 10*np.ones(N)
        self.pressure = np.ones(N)
        self.potential = -64*np.ones(N)

        self.masscharge = 246.55/0.2*np.ones(N) #c=masscharge*mass/Volume
        self.tension = 10*np.ones(N)
        self.vcrit = -64*np.ones(N)
        self.vexpansion = 0.05*np.ones(N)
        self.Acrit = 150*np.ones(N)

    def copy(self, system):

        self.num_cells = system.num_cells
        self.positions = system.positions.copy()
        # print("INSIDE COPY")
        # print(self.positions)
        # print(system.positions)
        # print("LEFT COPY")


        # These are the ion concentrations inside the cell
        self.ion_Cl = system.ion_Cl.copy()
        self.ion_K = system.ion_K.copy()
        self.ion_Na = system.ion_Na.copy()

        # prefactors in ion concentration. dry mass fraction ratio of channels vs atpase time alpha
        self.chanel_Cl = system.chanel_Cl.copy()
        self.chanel_K = system.chanel_K.copy()
        self.chanel_Na = system.chanel_Na.copy()

        self.mass = system.mass.copy()
        self.height = system.height.copy()
        self.pressure = system.pressure.copy()
        self.potential = system.potential.copy()

        self.masscharge = system.masscharge.copy()
        self.tension = system.tension.copy()
        self.vcrit = system.vcrit.copy()
        self.vexpansion = system.vexpansion.copy()
        self.Acrit = system.Acrit.copy()
