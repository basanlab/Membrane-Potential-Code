import numpy as np
import constants as constants
import vmpg.output as output

def fix_ions(system):
    #hardwireing concentrations:
    system.ion_Na = 2 * constants.MEDIUM_ION_Na * np.ones(system.num_cells)
    system.ion_K  = 2 * constants.MEDIUM_ION_K * np.ones(system.num_cells)
    system.ion_Cl = 2 * constants.MEDIUM_ION_Cl * np.ones(system.num_cells)


def fix_pressure(system):
    system.pressure = 5*np.ones(system.num_cells)

def constantnumber(system,area):
    for idx in range(system.num_cells):
        system.ion_Na[idx]= constants.MEDIUM_ION_Na * (1 + 1. / area[idx])
        system.ion_K[idx]= constants.MEDIUM_ION_K * (1 + 1. / area[idx])
        system.ion_Cl[idx]= constants.MEDIUM_ION_Cl * (1 + 1. / area[idx])

