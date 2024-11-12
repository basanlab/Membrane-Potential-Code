#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Living Matter Group - ICS-2 Juelich Research Center <j.elgeti@fz-juelich.de>
#
# Distributed under terms of the MIT license.


from vmpg import output
from vmpg import boundaries
from vmpg.system import System
import json
import numpy as np
from constants import *

#creation parameters
#N=56 #, L=75 works fine 7x8 cells grid
#ebenso 12, [75./7*3, 75./2]
#simulation_box = [75, 75]
N=50
simulation_box = [75./7*3, 75./2]
simulation_box = [200, 200]
initname='init50.vtu'

system = System(N,simulation_box)
#hex lattice:
Area=simulation_box[0]*simulation_box[1]
dx=np.sqrt(Area/N/np.sin(60.*np.pi/180.))
dy=np.sin(60.*np.pi/180.)*dx
x0=0
x=0
y=0
l=0
step=0
for idx in range(N):
    system.positions[idx]=np.array([x,y])
    x+=dx
    if(x>simulation_box[0]):
        y+=dy
        x0+=dx/2*(-1)**l
        l+=1 #alternate sign
        x=x0

#randomize positions:
for idx in range(N):
    system.positions[idx] += 0*np.random.rand(2)


boundaries.periodicbc(system,simulation_box)



output.write_vtu("outputs.pvd",initname, system, simulation_box,step, step*DT)

