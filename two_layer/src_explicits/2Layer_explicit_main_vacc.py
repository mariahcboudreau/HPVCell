#!/users/m/c/mcboudre/anaconda3/bin/python3
# HPVCellSim Master Equations VACC

import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt
from numba import jit
import sys
import time



# WE will use the odeint routine
from scipy.integrate import odeint
# With a wrapper to facilitate 2d arrays
from odeintw import odeintw

# Master Equations




@jit(nopython=True)
def J(c, t, beta, gamma, delta, rho, theta):
    """
        Time derivative of the occupation numbers.

            * c is the state distribution (array like)
            * t is time (scalar)
            * beta is the basal to basal/basal division rate
            * gamma is the basal to basal/parabasal division rate
            * delta is the basal to parabasal/parabasal division rate
            * rho is the parabasal to parabasal/parabasal divivsion rate
            * theta is the parabasal cell shedding/death rate
            * NB: We use logistic growth for preys to limit the # of states
            * K will be the carrying capacity
        """
    dx = 0 * c
    # dx[b, p] is the entire master equation
    for b, p in np.ndindex(c.shape):
        if b < (c.shape[0] - 1): # basal to basal/basal output
            dx[b, p] -= beta * b * c[b, p]
        if p < c.shape[1] - 1:  # basal to parabasal/basal output
            dx[b, p] -= gamma * b * c[b, p]
        if p < c.shape[1] - 1:  # parabasal to parabasal/parabasal output
            dx[b, p] -= rho * p * c[b, p]
        if p < c.shape[1] - 2:  # basal to parabasal/parabasal output
            dx[b, p] -= delta * b * c[b, p]
        if p > 0: # parabasal cell shedding output
            dx[b, p] -= theta * p * c[b, p]
        if b > 1 :  # basal to basal/basal input
            dx[b, p] += beta * (b - 1) * c[b - 1, p]
        if b > 0 and p > 0:  # basal to parabasal/basal input
            dx[b, p] += gamma * b * c[b, p - 1]
        if p > 1:  # parabasal to parabasal/parabasal input
            dx[b, p] += rho * (p - 1) * c[b, p - 1]
        if b < (c.shape[0] - 1) and p > 1:  # basal to parabasal/parabasal input
            dx[b, p] += delta * (b + 1) * c[b + 1, p - 2]
        if p < c.shape[1] - 1: # parabasal cell shedding input
            dx[b, p] += theta * (p + 1) * c[b, p + 1]

    return dx

# m[0] first moment for b
# m[1] first moment for p
# m[2] second moment for b
# m[3] second moment for p
# m[4] interaction term bp


#MoM
# def MOM(m, t, beta, gamma, delta, rho, theta): #figure out the right way to do the average of separate variables
#     dx = 0*m
#     # First Moment for b
#     dx[0] = (beta - delta)*m[0]
#     # First Moment for p
#     dx[1] = (rho - theta) * m[1] + (2 * delta + gamma) * m[0]
#     # Second Moment for b
#     dx[2] = 2*(beta - delta) * m[2] + (beta + delta) * m[0]
#     # Second Moment for p
#     dx[3] = (theta + rho) * m[1] + (2 * rho - 2 * theta) * m[3] + (2 * gamma + 4 * delta) * m[4] + (gamma + 4 * delta) * m[0]
#     # Covariance for bp
#     dx[4] = beta * m[4] - theta * m[4] + rho * m[4] + gamma * m[2] - delta * m[4] + 2 * delta * m[2] - 2 * delta * m[0]
#     return dx

#  Time declaration

time_total = 200
t_vec = np.linspace(0,time_total, time_total)


# Initial conditions
nb_of_states_b = 50
nb_of_states_p = nb_of_states_b * 2

x_0_delta = np.zeros((nb_of_states_b, nb_of_states_p))
x_0_poisson = np.zeros((nb_of_states_b, nb_of_states_p))
x_0_delta[1,0] = 1 # This is a delta approximation for the initial conditions x_0[1,0] = 1

for i in range(len(x_0_poisson)):
    x_0_poisson[i,0] = scipy.stats.poisson.pmf(k = i, mu = 1)

# m_0_delta = np.zeros(5)
# m_0_poisson = np.zeros(5)

# m_0_delta[0] = 1
# m_0_delta[2] = 1 

# m_0_poisson[0] = 1
# m_0_poisson[2] = 2

x_delta = x_0_delta
x_poisson = x_0_poisson
# m_delta = m_0_delta
# m_poisson = m_0_poisson

cumu_extinct_delta = np.zeros((time_total))
cumu_extinct_delta_last = np.zeros((time_total))
cumu_extinct_poisson_last = np.zeros((time_total))
# extinct_mom_b_delta = np.zeros((time_total))
# extinct_mom_b_poisson = np.zeros((time_total))

# Parameters of the model
symm_div = 0.04
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
shed = 2.0*(0.0082) # similar to division but just a little bit different according to the plos epithelial strat paper

rho = symm_div
gamma = asymm_div
delta = symm_div
beta = symm_div + 0.01
theta = shed



#######
###
### Time for loop
###
######

G = lambda x, t: J(x, t, beta, gamma, delta, rho, theta)

# for ti in range(0, time_total):

# Integration
print("integrating 1 ")
x_path_delta = odeintw(G, x_delta, t_vec)


print('integrating 2')
x_path_poisson = odeintw(G, x_poisson, t_vec)

# M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
# m_path_delta = odeintw(M, m_delta, t_vec)

# M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
# m_path_poisson = odeintw(M, m_poisson, t_vec)



# # Plot
# for t in range(0,time_total, 9):
#     #plt.plot(range(nb_of_states_b), np.sum(x_path_delta[t], axis=1), marker="o", ls='--', label=fr"$t = {t}$ for Delta initial conditions")
#     plt.plot(range(nb_of_states_b), np.sum(x_path_poisson[t], axis=1), marker="o", ls='--', label=fr"$t = {t}$ for Poisson initial conditions")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of particles')
# plt.show()
for t in range(0,time_total, 50):
    plt.plot(range(nb_of_states_p), np.sum(x_path_delta[t], axis=0), marker="o", ls='--', label=fr"$t = {t}$ for Delta initial conditions")
    #plt.plot(range(nb_of_states_b), np.sum(x_path_poisson[t], axis=1), marker="o", ls='--', label=fr"$t = {t}$ for Poisson initial conditions")
plt.legend()
plt.ylabel('Occupation number')
plt.xlabel('Number of particles')
plt.show()
    ####
    ## Solving for extinction probabilities 
    ####

        # Explicit ME

for ti in range(0, time_total):
    print("Calculating")
    cumu_extinct_delta[ti] = np.sum(x_path_delta[ti][0][:])

    cumu_extinct_delta_last[ti] = np.sum(x_path_delta[ti][:][0])
    cumu_extinct_poisson_last[ti] = np.sum(x_path_poisson[ti][:][0])
    
    # # # MOM
    # if m_path_delta[1][2] != 0:
    #     extinct_mom_b_delta[ti] = (1-((m_path_delta[ti][0]**2)/(m_path_delta[ti][2])))
    # else: 
    #     extinct_mom_b_delta[ti] = 0

    # if ( ((m_path_poisson[ti][2] - m_path_poisson[ti][0]) != 0) & (ti > 9) ):
    #     extinct_mom_b_poisson[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][2] - m_path_poisson[ti][0]) 
    # else: 
    #     extinct_mom_b_poisson[ti] = 0

    # print('reassignment')
    # x_delta = x_path_delta[1]
    # x_poisson = x_path_poisson[1]
    # m_delta = m_path_delta[1]
    # m_poisson = m_path_poisson[1]

with open('cumu_extinct_2layer_main_delta_last_12-11.npy', 'wb') as handle:
    np.save(handle, cumu_extinct_delta_last)

with open('cumu_extinct_2layer_main_poisson_last_12-11.npy', 'wb') as handle:
    np.save(handle, cumu_extinct_poisson_last)

plt.plot(t_vec, cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta Approx)')
# plt.plot(t_vec, cumu_extinct_poisson, label = 'Cumulative probability of extinction - Explicit basals (Poisson Approx)')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('2 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()







# print('saving')
# sys.stdout.flush()
# with open('/gpfs1/home/m/c/mcboudre/scratch/hpv_modeling/cumu_extinct_2layer_main_delta.npy', 'wb') as handle:
#     np.save(handle, cumu_extinct_delta)
# with open('/gpfs1/home/m/c/mcboudre/scratch/hpv_modeling/cumu_extinct_2layer_main_poisson.npy', 'wb') as handle:
#     np.save(handle, cumu_extinct_poisson)

# with open('/gpfs1/home/m/c/mcboudre/scratch/hpv_modeling/extinct_mom_b_2layer_main_delta_loop.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_delta)
# with open('/gpfs1/home/m/c/mcboudre/scratch/hpv_modeling/extinct_mom_b_2layer_main_poisson_loop.npy', 'wb') as handle:    
#     np.save(handle, extinct_mom_b_poisson)



print('done')
# sys.stdout.flush()