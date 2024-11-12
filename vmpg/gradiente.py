
import numpy as np

from vmpg import voronoi


def calcgradient(system,simulation_box):
    forces = np.zeros_like(system.positions)
    N = len(system.positions)
    pos=system.positions
    eps=0.001
    E0=voronoi.calculate_voronoi_energy(system.positions,system,simulation_box)
    #looping over particles idx
    for idx in range(N):
        #manually calculating gradient brute force
        #force = -dE/dx, minus sign order of terms
        pos[idx][0]+=eps
        forces[idx][0] = (E0-voronoi.calculate_voronoi_energy(pos,system,simulation_box))/eps
        pos[idx][0]-=eps
        pos[idx][1]+=eps
        forces[idx][1] = (E0-voronoi.calculate_voronoi_energy(pos,system,simulation_box))/eps
        pos[idx][1]-=eps

    return forces
