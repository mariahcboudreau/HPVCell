######
## HPVCellSim Master Equations MOM 5 layer with Dead Cells
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
#       m[2]  - first moment for i
#       m[3]  - first moment for s
#       m[4]  - first moment for d
#       m[5]  - second moment for b
#       m[6]  - second moment for p
#       m[7]  - second moment for i
#       m[8]  - second moment for s
#       m[9]  - second moment for d
#       m[10] - interation term bp
#       m[11] - interation term bi
#       m[12] - interation term bs
#       m[13] - interation term bd
#       m[14] - interation term pi
#       m[15] - interation term ps
#       m[16] - interation term pd
#       m[17] - interation term is
#       m[18] - interation term id
#       m[19] - interation term sd
#       t     - time
#       beta  - Rate of bbb division
#       gamma - Rate of bbp division
#       delta - Rate of bpp division
#       rho   - Rate of ppp division
#       alpha - Rate of pi movement
#       sigma - Rate of is movement
#       theta - Rate of shed 
#   Outputs: 
#       Updated m values 
###################################

def MOM(m, t, beta, gamma, delta, rho, alpha, sigma, theta): 
    dx = 0*m
    # First Moment for b
    dx[0] = (beta - delta) * m[0]

    # First Moment for p
    dx[1] = (rho - alpha) * m[1] + (2 * delta + gamma) * m[0]

    # First Moment for i
    dx[2] = alpha * m[1] - sigma * m[2]

    # First Moment for s
    dx[3] = sigma * m[2] - theta * m[3]

    # First Moment for d
    dx[4] = theta * m[3]

    # Second Moment for b
    dx[5] = 2 * (beta - delta) * m[5] + (beta + delta) * m[0]

    # Second Moment for p
    dx[6] = 2 * (rho - alpha) * m[6] + (rho + alpha) * m[1] + 2 * (gamma + 2 * delta) * m[10] + (gamma + 4 * delta) * m[0]

    # Second Moment for i
    dx[7] = -2 * sigma * m[7] + sigma * m[2] + alpha * m[1] + 2 * alpha * m[14]

    # Second Moment for s
    dx[8] = -2 * theta * m[8] + theta * m[3] + sigma * m[2] + 2 * sigma * m[17]

    # Second Moment for d
    dx[9] = 2 * theta * m[19] + theta * m[3]

    # Covariance for bp
    dx[10] = (beta + rho - delta - alpha) * m[10] + (gamma + 2 * delta) * m[5] - 2 * delta * m[0]

    # Covariance for bi
    dx[11] = (beta - delta - sigma) * m[11] + alpha * m[10]

    # Covariance for bs
    dx[12] = (beta - delta - theta) * m[12] + sigma * m[11]

    # Covariance for bd
    dx[13] = (beta - delta) * m[13] + theta * m[12]

    # Covariance for pi
    dx[14] = (rho - alpha - sigma) * m[14] + (gamma + 2 * delta) * m[11] + alpha * m[6] - alpha * m[1]

    # Covariance for ps
    dx[15] = (rho - theta - alpha) * m[15] + (gamma + 2 * delta) * m[12] + sigma * m[14]

    # Covariance for pd
    dx[16] = (rho - alpha) * m[16] + (gamma + 2 * delta) * m[13] + theta * m[15]

    # Covariance for is
    dx[17] = -1*(theta + sigma) * m[17] + alpha * m[15] + alpha * m[7] - alpha * m[2]

    # Covariance for id
    dx[18] = theta * m[17] - sigma * m[18] + alpha * m[16] 

    # Covariance for sd
    dx[19] = theta * m[8] - theta * m[19] - theta * m[3] + sigma * m[18]

    return dx


# Max Time
t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions

m_0_delta_geo = np.zeros(20)
m_0_delta_geo[0] = 1
m_0_delta_geo[5] = 1

# m_0_poisson = np.zeros(9)
# m_0_poisson[0] = 1
# m_0_poisson[3] = 2


# Updated divsion parameters
R_b = 0.03          # Division rate of basal - Murall citation
R_p = 0.39          # Division rate of parabasal - Murall citation
symm_div = 0.08     # Symmetric division rate - Clayton
asymm_div = 0.84    # Asymmetric division rate - Clayton

beta = R_b * symm_div  + 0.001      # bbb 
gamma = R_b * asymm_div             # bbp
delta = R_b * symm_div              # bpp
rho = R_p * symm_div                # ppp
alpha = 0.4                         # transition pi - Murall
sigma = 0.4                         # transition is - Murall
theta = 0.67                        # Shed - Murall


################################# INTEGRATION ###################################

M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, alpha, sigma, theta)
# m_path_poisson = odeintw(M, m_0_poisson, t_vec)
m_path_delta_geo = odeintw(M, m_0_delta_geo, t_vec)



################################# EXTINCTION DERIVATION ###################################

# Creating extinction variables

# extinct_mom_b_delta = np.zeros((len(t_vec)))
# extinct_mom_b_poisson = np.zeros((len(t_vec)))
# extinct_mom_p_poisson = np.zeros((len(t_vec)))
extinct_mom_b_geometric = np.zeros((len(t_vec)))
extinct_mom_p_geometric = np.zeros((len(t_vec)))
# extinct_mom_b_expon = np.zeros((len(t_vec)))
# extinct_mom_b_power = np.zeros((len(t_vec)))


for ti in range(0, len(t_vec)):

    # Poisson and geometric approximations

    # extinct_mom_b_poisson[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][3] - m_path_poisson[ti][0])  + scipy.stats.poisson.pmf(k = 0, mu = ((m_path_poisson[ti][3])/(m_path_poisson[ti][0])-1))
    # extinct_mom_p_poisson[ti] = 1 - (m_path_poisson[ti][1]**2)/(m_path_poisson[ti][4] - m_path_poisson[ti][1])  + scipy.stats.poisson.pmf(k = 0, mu = ((m_path_poisson[ti][4])/(m_path_poisson[ti][1])-1))
    extinct_mom_b_geometric[ti] = 1 - (2*m_path_delta_geo[ti][0]**2) / (m_path_delta_geo[ti][5] + m_path_delta_geo[ti][0] )  + scipy.stats.geom.pmf(k = 0, p = (2/(m_path_delta_geo[ti][0]*m_path_delta_geo[ti][5] + 1)))
    # extinct_mom_p_geometric[ti] = 1 - (2*m_path_delta_geo[ti][1]**2) / (m_path_delta_geo[ti][3] + m_path_delta_geo[ti][1] )  + scipy.stats.geom.pmf(k = 0, p = (2/(m_path_delta_geo[ti][1]*m_path_delta_geo[ti][4] + 1)))

############################## SAVING FILES #######################################


# with open('data/extinction_mom_b_2layerdead_poisson_'+date+'_time500.npy', 'wb') as f:
#     np.save(f, extinct_mom_b_poisson)
with open('four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'wb') as f:
    np.save(f, extinct_mom_b_geometric)

with open('four_layer/data/shed_first_moments_delta_4layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,4])

with open('four_layer/data/basal_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,0])

with open('four_layer/data/para_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,1])

with open('four_layer/data/intermed_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,2])

with open('four_layer/data/super_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,3]) 

with open('four_layer/data/para_second_moment_geom_4layerwithdead_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,6])



# with open('extinct_mom_b_2layer_main_geometric_1_8_500.npy', 'wb') as handle:    
#     np.save(handle, extinct_mom_b_geometric)
# with open('extinct_mom_b_2layer_main_delta_11_20.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_delta_geo)

# with open('extinct_mom_b_2layer_main_geometric_11_20.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_geometric)



############################## GET FUNCTIONS #######################################

# def get_m_path_poisson():
#     return m_path_poisson

# def get_m_path_poisson():
#     return m_path_poisson
