import numpy as np
from scipy.optimize import broyden1
import constants as constants
from numba import jit

def calculate_ions_steady_state(system, area, volume):
    """ system is the system, area is the surface area of the cells  """
    #looping over cells, updating ions and potential. To Optimize: Vectorize
    for idx in range(system.num_cells):
        system.ion_Cl[idx], system.ion_K[idx], system.ion_Na[idx],system.potential[idx] = \
            solve_ions_steady_state_cell_dynamic(volume[idx], area[idx], system.masscharge[idx]*system.mass[idx],
                                                 system.chanel_K[idx], system.chanel_Na[idx],constants.SATURATION,
                                                 system.ion_Cl[idx], system.ion_K[idx], system.ion_Na[idx])
#        system.ion_Cl[idx], system.ion_K[idx], system.ion_Na[idx],system.potential[idx] = \
#            solve_ions_steady_state_cell(volume[idx], area[idx], system.mass[idx],
#                                                 system.chanel_K[idx], system.chanel_Na[idx],
#                                                 system.ion_Cl[idx], system.ion_K[idx], system.ion_Na[idx])




@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def solve_ions_steady_state_cell_dynamic(volume, area, massc, channel_K, channel_Na,sat, ion_Cl, ion_K, ion_Na):
    mc=massc/volume
    #getting a guess from qubic fit:
    #Hallo anpassen
    #base:
    u=-26.181786589278243+-0.26539874891505294*mc+0.0005130875533623571*mc**2+-4.5781042575789454e-07*mc**3
    ion_Na=4.1624955363354195+0.006452532120880027*mc+-1.2978630531635247e-05*mc**2+1.1808248896673846e-08*mc**3

    
    #val1 BchK=BchK+BchK*1
    #u=-32.27637692499921+-0.3207628441671858*mc+0.0007096954591550898*mc**2+-6.944516351622768e-07*mc**3
    #ion_Na=4.313562205457082+0.007693372232416019*mc+-1.7653440679028816e-05*mc**2+1.7547811628083072e-08*mc**3


    #gram0p2 model='gram0p5'BchK=BchK+BchNa*0.2 BchNa=BchNa+BchNa*0.2
    #u=-24.660525432750003+-0.2509995740008313*mc+0.00046137850511147654*mc**2+-3.940523841397053e-07*mc**3
    #ion_Na=4.4024166919151355+0.006648147363821969*mc+-1.2580864997765941e-05*mc**2+1.0900879912109756e-08*mc**3


    #than optain Cl and K following those:
    ion_Cl=constants.MEDIUM_ION_Cl*np.exp(u*constants.FARADAY/constants.RG/constants.TEMP)
    ion_K=ion_Cl+mc-ion_Na #use K here, that is always high enough it does not go negative


    def iterate(dt,K,Na):
        Cl= Na+K-mc   #note: Na+K>c M!!!
        mu = np.log(Cl/constants.MEDIUM_ION_Cl)
        emu = np.exp(mu)
        delNa = channel_Na * mu * (constants.MEDIUM_ION_Na-Na*emu)/(1-emu)+1*(Na**3/(sat**3+Na**3))#multiplied by MM kinetic to avoid negative insede
        delK = channel_K * mu * (constants.MEDIUM_ION_K-K*emu)/(1-emu)-2/3*(Na**3/(sat**3+Na**3))#multiplied by MM kinetic, na in limited
        Na -= dt*delNa
        K -= dt*delK
        K=max(K,mc-Na+5) #avoids negative CL, with +5 limits to som 90mv or so
        #print(delK,delNa)
        return [Cl,K,Na,mu,delNa**2+delK**2]
    delta=1
    zeitschritt=0.0005*sat**3+0.01
    #print('--------------------------------------------')
    #print(volume, area, massc, channel_K, channel_Na,sat, ion_Cl, ion_K, ion_Na)
    #print(ion_Cl, ion_K, ion_Na, u,iters)
    while delta>1e-10 :  #e-14 is still stable and more accurate. Lets keep 1e-10 for now, which fits with the estimator... 
        #print('concentrationsb: ',ion_Cl, ion_K, ion_Na,delta,zeitschritt)
        ion_Cl, ion_K, ion_Na, mu,delta = iterate(zeitschritt,ion_K, ion_Na)
        #if(delta>pdelta):
        #    zeitschritt*=0.5
        #print('concentrations: ',ion_Cl, ion_K, ion_Na,mu,delta,zeitschritt)
    u = mu/constants.FARADAY*constants.RG*constants.TEMP
    return ion_Cl, ion_K, ion_Na, u





#alternatively, the INVERSE Function is analytic. 
#i.e. calculating the channel_Na, masscharge, etc from the ion concentrations. 
import numpy as np
def channelofion(Clj,Kj,Naj,sat):
    muj=np.log(Clj/constants.MEDIUM_ION_Cl)
    uj = muj/FARADAY*RG*TEMP
    emuj =  np.exp(muj)  
    #charge of drymassdensity given by charge balance. no Constant => inside Cl+bicarbonate
    drymasscharge=Naj+Kj-Clj
    #chanel factors given by potential and concentrations
    channel_Na = -1. * (1-emuj)/ muj / (constants.MEDIUM_ION_Na-Naj*emuj)*(Naj**3/(sat**3+Naj**3))
    channel_K = 2./3 * (1-emuj)/ muj / (constants.MEDIUM_ION_K-Kj*emuj)*(Naj**3/(sat**3+Naj**3))
    return(drymasscharge,channel_K,channel_Na,uj)
