MEDIUM_ION_Cl = 160.0  #in miliMolar   
MEDIUM_ION_K = 5.0   
MEDIUM_ION_Na = 155.0  
FARADAY = 0.096485 # 96485/1000./1000.0 #we are using mili molar, and also mV  
RG = 0.008314 #8.314/1000 #we are using mili molar  
TEMP = 300  
DT = 0.05  
GAMMA_BACKGROUND = 1000 #sets the pressure relaxation time.  
CELL_OSMOLITES = 80 #addtional osmolites in the cell, compared to medium. Adds to pressure  
SATURATION = 10 #Satturation constant of Na-K-Atpase 
simulation_box = [200.0, 200.0]  
dump_every = 20  
num_steps = 1000 #run MUCH longer to reproduce paper data  
output_name = "cells"  
inputfile = "init50.vtu" #will restart from that file!  
curstep = 0  
