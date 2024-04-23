######
## HPVCellSim Master Equations MOM 3 layer with Dead Cells
######


################################# IMPORTS ###################################
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import sys
import scipy
from datetime import datetime
date = datetime.today().strftime('%m-%d-%Y')

# WE will use the odeint routine
from scipy.integrate import odeint
# With a wrapper to facilitate 2d arrays
from odeintw import odeintw

##################################
# MOM Function
#   Description: Defines the governing functions of the cell division process for HPV
#  
#   Inputs: 
#       m[0]  - first moment for b
#       m[1]  - first moment for p
#       m[2]  - first moment for d
#       m[3]  - second moment for b
#       m[4]  - second moment for p
#       m[5]  - second moment for d
#       m[6]  - interaction term bp
#       m[7]  - interaction term pd
#       m[8]  - interaction term bd
#       t     - time
#       beta  - Rate of bbb division
#       gamma - Rate of bbp division
#       delta - Rate of bpp division
#       rho   - Rate of ppp division
#       theta - Rate of shed 
#   Outputs: 
#       Updated m values 
###################################

def MOM(m, t, beta, gamma, delta, rho, theta): 
    dx = 0*m
    # First Moment for b
    dx[0] = (beta - delta)* m[0]
    # First Moment for p
    dx[1] = (rho - theta) * m[1] + (2 * delta + gamma) * m[0]
    # First Moment for d
    dx[2] = theta * m[1]
    # Second Moment for b
    dx[3] = 2*(beta - delta) * m[3] + (beta + delta) * m[0]
    # Second Moment for p
    dx[4] = (theta + rho) * m[1] + (2 * rho - 2 * theta) * m[4] + (2 * gamma + 4 * delta) * m[6] + (gamma + 4 * delta) * m[0]
    # Second Moment for d
    dx[5] = 2*theta * m[7] + theta * m[1]
    # Covariance for bp
    dx[6] = beta * m[6] - theta * m[6] + rho * m[6] + gamma * m[3] - delta * m[6] + 2 * delta * m[3] - 2 * delta * m[0]
    # Covariance for pd
    dx[7] = (rho - theta) * m[7] + theta * m[4] - theta * m[1] + (gamma + 2*delta) * m[8]
    # Covariance for bd
    dx[8] = beta * m[0] + theta * m[6] - delta * m[8]
    return dx


# Max Time
t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions

m_0_delta_geo = np.zeros(9)
m_0_delta_geo[0] = 1
m_0_delta_geo[3] = 1

m_0_poisson = np.zeros(9)
m_0_poisson[0] = 1
m_0_poisson[3] = 2


# Updated divsion parameters
R_b = 0.03          # Division rate of basal - Murall citation
R_p = 0.39          # Division rate of parabasal - Murall citation
symm_div = 0.08     # Symmetric division rate - Clayton
asymm_div = 0.84    # Asymmetric division rate - Clayton

beta = R_b * symm_div  + 0.001      # bbb 
gamma = R_b * asymm_div             # bbp
delta = R_b * symm_div              # bpp
rho = R_p * symm_div                # ppp
theta = 0.67                        # Shed - Murall


################################# INTEGRATION ###################################

M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
m_path_poisson = odeintw(M, m_0_poisson, t_vec)
m_path_delta_geo = odeintw(M, m_0_delta_geo, t_vec)



################################# EXTINCTION DERIVATION ###################################

# Creating extinction variables

# extinct_mom_b_delta = np.zeros((len(t_vec)))
extinct_mom_b_poisson = np.zeros((len(t_vec)))
extinct_mom_p_poisson = np.zeros((len(t_vec)))
extinct_mom_b_geometric = np.zeros((len(t_vec)))
extinct_mom_p_geometric = np.zeros((len(t_vec)))
# extinct_mom_b_expon = np.zeros((len(t_vec)))
# extinct_mom_b_power = np.zeros((len(t_vec)))


for ti in range(0, len(t_vec)):

    # extinct_mom_b_delta[ti] = (1-((m_path_delta_geo[ti][0]**2)/(m_path_delta_geo[ti][2])))
    
    # Poisson and geometric approximations

    extinct_mom_b_poisson[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][3] - m_path_poisson[ti][0])  + scipy.stats.poisson.pmf(k = 0, mu = ((m_path_poisson[ti][3])/(m_path_poisson[ti][0])-1))
    extinct_mom_p_poisson[ti] = 1 - (m_path_poisson[ti][1]**2)/(m_path_poisson[ti][4] - m_path_poisson[ti][1])  + scipy.stats.poisson.pmf(k = 0, mu = ((m_path_poisson[ti][4])/(m_path_poisson[ti][1])-1))
    extinct_mom_b_geometric[ti] = 1 - (2*m_path_delta_geo[ti][0]**2) / (m_path_delta_geo[ti][3] + m_path_delta_geo[ti][0] )  + scipy.stats.geom.pmf(k = 0, p = (2/(m_path_delta_geo[ti][0]*m_path_delta_geo[ti][3] + 1)))
    extinct_mom_p_geometric[ti] = 1 - (2*m_path_delta_geo[ti][1]**2) / (m_path_delta_geo[ti][3] + m_path_delta_geo[ti][1] )  + scipy.stats.geom.pmf(k = 0, p = (2/(m_path_delta_geo[ti][1]*m_path_delta_geo[ti][4] + 1)))
    # shed_moments[0][ti] = 1000*theta*m_path_delta_geo[ti][1]
    # shed_moments[1][ti] = 1000**2*theta*m_path_delta_geo[ti][3]
    # # extinct_mom_b_expon[ti] = 1 - (2*m_path_poisson[ti][0]**2) / (m_path_poisson[ti][2]) #+ scipy.stats.expon.pdf(x = 0, scale = 1/((2*m_path_poisson[ti][0])/(m_path_poisson[ti][2])))  
    # extinct_mom_b_power[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][2]) + scipy.stats.powerlaw.pdf(x=0, a = ((3*m_path_poisson[ti][2] - 2*m_path_poisson[ti][0])/(m_path_poisson[ti][2]-m_path_poisson[ti][0])))


############################## SAVING FILES #######################################

import os
os.chdir('two_layer/src')
# with open('data/extinction_mom_b_2layerdead_poisson_'+date+'_time500.npy', 'wb') as f:
#     np.save(f, extinct_mom_b_poisson)
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time500.npy', 'wb') as f:
    np.save(f, extinct_mom_b_geometric)
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/shed_first_moments_delta_2layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,2])
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/basal_first_moment_geom_2layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,0])
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/para_first_moment_geom_2layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,1])
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/para_second_moment_geom_2layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,4])
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/basal_second_moment_geom_2layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,3])



# with open('extinct_mom_b_2layer_main_geometric_1_8_500.npy', 'wb') as handle:    
#     np.save(handle, extinct_mom_b_geometric)
# with open('extinct_mom_b_2layer_main_delta_11_20.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_delta_geo)

# with open('extinct_mom_b_2layer_main_geometric_11_20.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_geometric)



############################## GET FUNCTIONS #######################################

def get_m_path_poisson():
    return m_path_poisson

def get_m_path_poisson():
    return m_path_poisson

