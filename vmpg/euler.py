# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:50:12 2019

@author: tobias
"""
import constants as constants


# Here we suppose we are in the overdamped regime. That is, the 
# inercia play no role in the system
def advance_positions(system, forces):
    system.positions += forces / (constants.GAMMA_BACKGROUND) * constants.DT
    #Frage: warum mas, np newaxis?
