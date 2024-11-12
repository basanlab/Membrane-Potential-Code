#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Living Matter Group - ICS-2 Juelich Research Center <j.elgeti@fz-juelich.de>
#
# Distributed under terms of the MIT license.

import numpy as np
from constants import *
from vmpg import boundaries

def calculate_osmotic_pressure(system):
    system.pressure = RG*TEMP*(
        (system.ion_Na - MEDIUM_ION_Na) +
        (system.ion_K  - MEDIUM_ION_K ) +
        (system.ion_Cl - MEDIUM_ION_Cl)
        + CELL_OSMOLITES
    )

