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
t_length = 200
t_steps = 200
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
beta = symm_div + 0.01
theta = shed



# Integration
M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
m_path_poisson = odeintw(M, m_0_poisson, t_vec)
m_path_delta = odeintw(M, m_0_delta, t_vec)

print("done integrating")

extinct_mom_b_delta = np.zeros((len(t_vec)))
extinct_mom_b_poisson = np.zeros((len(t_vec)))
   # # MOM
for ti in range(0, len(t_vec)):

    if m_path_delta[ti][2] != 0:
        extinct_mom_b_delta[ti] = (1-((m_path_delta[ti][0]**2)/(m_path_delta[ti][2])))
    else: 
        extinct_mom_b_delta[ti] = 0

    # if ( ((m_path_poisson[ti][2] - m_path_poisson[ti][0]) != 0) & (ti > 9) ):
    extinct_mom_b_poisson[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][2] - m_path_poisson[ti][0]) 
# else: 
    # extinct_mom_b_poisson[ti] = 0

plt.plot(t_vec, extinct_mom_b_delta, label = 'Cumulative probability of extinction - MOM basals (Delta Approx)')
plt.plot(t_vec, extinct_mom_b_poisson, label = 'Cumulative probability of extinction - MOM basals (Poisson Approx)')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('2 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()

# plt.plot(t_vec, extinct_mom_b_poisson, label = 'Cumulative probability of extinction - MOM basals (Poisson Approx)')
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# plt.title('2 Layer system')
# #plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()

# print('almost done pickling')
# with open('extinct_mom_b_2layer_main_poisson_11_17.npy', 'wb') as handle:    
#     np.save(handle, extinct_mom_b_poisson)
# with open('extinct_mom_b_2layer_main_delta_11_17.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_delta)