import numpy as np


def periodicbc(system,simulation_box):
    for idx in range(system.num_cells):
        #??? write this more elegantly
        system.positions[idx][0]-=simulation_box[0]*np.floor(system.positions[idx][0]/simulation_box[0])
        system.positions[idx][1]-=simulation_box[1]*np.floor(system.positions[idx][1]/simulation_box[1])

def minimalr(r,simulation_box):
    r[0]-=simulation_box[0]*np.round(r[0]/simulation_box[0])
    r[1]-=simulation_box[1]*np.round(r[1]/simulation_box[1])
    return r
