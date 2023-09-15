# HPVCellSim Master Equations - 4 Layer

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import sys



# WE will use the odeint routine
from scipy.integrate import odeint
# With a wrapper to facilitate 2d arrays
from odeintw import odeintw


@jit(nopython=True)


#MoM
def MOM(m, t, beta, gamma, delta, rho, alpha, sigma, theta): #figure out the right way to do the average of separate variables
    dx = 0*m
    # First Moment for b
    dx[0] = (beta - delta) * m[0]

    # First Moment for p
    dx[1] = (rho - alpha) * m[1] + (2 * delta + gamma) * m[0]

    # First Moment for i
    dx[2] = -sigma * m[2]

    # First Moment for s
    dx[3] = -theta * m[3]

    # Second Moment for b
    dx[4] = 2 * (beta - delta) * m[4] + (beta + delta) * m[0]

    # Second Moment for p
    dx[5] = 2 * (rho - alpha) * m[5] + (rho + alpha) * m[1] + 2 * (gamma + 2 * delta) * m[8] + (gamma + 4 * delta) * m[0]

    # Second Moment for i
    dx[6] = -2 * sigma * m[6] + sigma * m[2] 

    # Second Moment for s
    dx[7] = -2 * theta * m[7] + theta * m[3]

    # Covariance for bp
    dx[8] = (beta + rho - delta - alpha) * m[8] + (gamma + 2 * delta) * m[4] - 2 * delta * m[0]

    # Covariance for bi
    dx[9] = (beta - delta - sigma) * m[9]

    # Covariance for bs
    dx[10] = (beta - delta - theta) * m[10]

    # Covariance for pi
    dx[11] = (rho - alpha - sigma) * m[11] + (gamma + 2 * delta) * m[9]

    # Covariance for ps
    dx[12] = (rho - theta - alpha) * m[12] + (gamma + 2 * delta) * m[10]

    # Covariance for is
    dx[13] = -1*(theta + sigma) * m[13]

    return dx


# Time of observations
t_length = 500
t_steps = 50
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions

m_0 = np.zeros(14)
m_0[0] = 1
m_0[4] = 1


# Parameters of the model
symm_div = 0.04
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
diff = 2.0*(0.0082)# similar to division but just a little bit different according to the plos epithelial strat paper

rho = symm_div
gamma = asymm_div
delta = symm_div
beta = symm_div - 0.01
theta = 0.01 # plost epithelial chart
alpha = diff
sigma = 0.18 # plos epithelial chart 




# Integration
M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, alpha, sigma, theta)
m_path = odeintw(M, m_0, t_vec)

 

print('almost done saving')
with open('m_path_4layer_main_1.npy', 'wb') as handle:    
    np.save(handle, m_path)
