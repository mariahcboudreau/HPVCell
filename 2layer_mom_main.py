# HPVCellSim Master Equations MOM 4 layer

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import sys



# WE will use the odeint routine
from scipy.integrate import odeint
# With a wrapper to facilitate 2d arrays
from odeintw import odeintw

# MOM

# m[0] first moment for b
# m[1] first moment for p
# m[2] second moment for b
# m[3] second moment for p
# m[4] interaction term bp


#MoM
def MOM(m, t, beta, gamma, delta, rho, theta): #figure out the right way to do the average of separate variables
    dx = 0*m
    # First Moment for b
    dx[0] = (beta - delta)* m[0]
    # First Moment for p
    dx[1] = (rho - theta) * m[1] + (2 * delta + gamma) * m[0]
    # Second Moment for b
    dx[2] = 2*(beta - delta) * m[2] + (beta + delta) * m[0]
    # Second Moment for p
    dx[3] = (theta + rho) * m[1] + (2 * rho - 2 * theta) * m[3] + (2 * gamma + 4 * delta) * m[4] + (gamma + 4 * delta) * m[0]
    # Covariance for bp
    dx[4] = beta * m[4] - theta * m[4] + rho * m[4] + gamma * m[2] - delta * m[4] + 2 * delta * m[2] - 2 * delta * m[0]
    return dx


# Time of observations
t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions

m_0_delta = np.zeros(5)
m_0_delta[0] = 1
m_0_delta[2] = 1

m_0_poisson = np.zeros(5)
m_0_poisson[0] = 1
m_0_poisson[2] = 2

# Parameters of the model
symm_div = 0.04
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
shed = 2.0*(0.0082)*1.99 # similar to division but just a little bit different according to the plos epithelial strat paper

rho = symm_div
gamma = asymm_div
delta = symm_div
beta = symm_div - 0.01
theta = shed



# Integration
M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
m_path_poisson = odeintw(M, m_0_poisson, t_vec)
m_path_delta = odeintw(M, m_0_delta, t_vec)





print('almost done pickling')
with open('extinct_mom_b_2layer_main_poisson.npy', 'wb') as handle:    
    np.save(handle, m_path_poisson)
with open('extinct_mom_b2layer_main_delta.npy', 'wb') as handle:
    np.save(handle, m_path_delta)