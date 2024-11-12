import numpy as np
from vmpg import system
import constants
from numpy.random import default_rng


def divideanddiebymassPot(areas,oldsystem,simulation_box,rng,tote):
    #finding new number of cells
    #rng = default_rng()
    oldn=oldsystem.num_cells
    newn=oldn
    dead = np.zeros(oldn)
    masses=oldsystem.mass
    pots=oldsystem.potential
    potC=oldsystem.vcrit
    for cell in range(oldn):
        #if(np.random.uniform()/constants.DT<(0.0057762265+ 0.0057762265*(pots[cell]<1.*potC[cell]))):
        if(rng.random()/constants.DT<(0.0057762265+ 0.0057762265*(pots[cell]<1.*potC[cell]))):
            #log(2)/100/24=0.000288811325233311 ln(2)/5/24= 0.0057762265
            tote+=1
            newn-=1
            dead[cell]=1
        if(masses[cell]>oldsystem.Acrit[cell] and not dead[cell]): #Acrit=critical mass
            newn+=1
    #print(pots,dead)
    newsystem = system.System(newn,simulation_box)
    count=0
    for cell in range(oldn):
        if(dead[cell]):
            continue
        copycell(oldsystem,newsystem,cell,count)
        count+=1
        if(masses[cell]>oldsystem.Acrit[cell] and not dead[cell]):
            copycell(oldsystem,newsystem,cell,count)
            theta=np.random.uniform(0,2*np.pi)
            newsystem.positions[count]+=0.1*np.array([np.cos(theta),np.sin(theta)]) #simple division in epsilonenv!
            newsystem.mass[count]*=0.5  #mass gets spread 50-50
            newsystem.mass[count-1]*=0.5
            #Acrit follows adder modell:
            newsystem.Acrit[count]=newsystem.mass[count]+rng.normal(loc=200,scale=50)  #Adder model, sigma=mean/2
            newsystem.Acrit[count-1]=newsystem.mass[count-1]+rng.normal(loc=200,scale=50)  #Adder model, sigma=mean/2
            count+=1
    return newsystem,tote


def copycell(fromsystem,tosystem,fromidx,toidx):
    tosystem.positions[toidx]=fromsystem.positions[fromidx]
    tosystem.ion_Cl[toidx]=fromsystem.ion_Cl[fromidx]
    tosystem.ion_K[toidx]=fromsystem.ion_K[fromidx]
    tosystem.ion_Na[toidx]=fromsystem.ion_Na[fromidx]
    tosystem.chanel_Cl[toidx]=fromsystem.chanel_Cl[fromidx]
    tosystem.chanel_K[toidx]=fromsystem.chanel_K[fromidx]
    tosystem.chanel_Na[toidx]=fromsystem.chanel_Na[fromidx]
    tosystem.mass[toidx]=fromsystem.mass[fromidx]
    tosystem.height[toidx]=fromsystem.height[fromidx]
    tosystem.pressure[toidx]=fromsystem.pressure[fromidx]
    tosystem.potential[toidx]=fromsystem.potential[fromidx]
    tosystem.masscharge[toidx]=fromsystem.masscharge[fromidx]
    tosystem.tension[toidx]=fromsystem.tension[fromidx]
    tosystem.vcrit[toidx]=fromsystem.vcrit[fromidx]
    tosystem.vexpansion[toidx]=fromsystem.vexpansion[fromidx]
    tosystem.vcrit[toidx]=fromsystem.vcrit[fromidx]
    tosystem.Acrit[toidx]=fromsystem.Acrit[fromidx]

    return 0

def grow(system):
    n=system.num_cells
    pots=system.potential
    for cell in range(n):
        system.mass[cell]*= 1 + (-0.0041258760+ (0.0041258760+0.0057762265)/(-61.1916849-system.vcrit[cell])*(system.potential[cell]-system.vcrit[cell])*(system.potential[cell]>system.vcrit[cell])) * constants.DT
        #note that membrane potential should be negative! thus minus minus
        #note that prefactor is rate/mV
        #Protein halflife about 7days=7*24h. Rate=(-ln(2)/7/24)=-0.0041258760
        #growth of DMD when empty: Doubles in 4days: Rate=(-ln(2)/4/24)=0.00722028313
        #sum for equation 0.0041258760+0.00722028313=0.0113461592\
        #changed growth starts at vcrit, and than linear to compensating apoptosis 0.0057762265 at our baseline (-61.1916849)
