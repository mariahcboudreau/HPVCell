#!/users/m/c/mcboudre/anaconda3/bin/python3
# HPVCellSim Master Equations VACC


import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt
from numba import jit
import sys




# We will use the odeint routine
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
def MOM(m, t, beta, gamma, delta, rho, theta): #figure out the right way to do the average of separate variables
    dx = 0*m
    # First Moment for b
    dx[0] = (beta - delta)*m[0]
    # First Moment for p
    dx[1] = (rho - theta) * m[1] + (2 * delta + gamma) * m[0]
    # Second Moment for b
    dx[2] = 2*(beta - delta) * m[2] + (beta + delta) * m[0]
    # Second Moment for p
    dx[3] = (theta + rho) * m[1] + (2 * rho - 2 * theta) * m[3] + (2 * gamma + 4 * delta) * m[4] + (gamma + 4 * delta) * m[0]
    # Covariance for bp
    dx[4] = beta * m[4] - theta * m[4] + rho * m[4] + gamma * m[2] - delta * m[4] + 2 * delta * m[2] - 2 * delta * m[0]
    return dx

#  Time declaration

time_total = 200
t_points = 20
t = np.linspace(0,time_total, t_points)


# Initial conditions
nb_of_states_b = 200
nb_of_states_p = nb_of_states_b * 2

x_0_poisson = np.zeros((nb_of_states_b, nb_of_states_p))

for i in range(len(x_0_poisson)):
    x_0_poisson[i,0] = scipy.stats.poisson.pmf(k = i, mu = 1)



x_poisson = x_0_poisson



cumu_extinct_poisson = np.zeros((time_total))



# Parameters of the model
symm_div = 0.04
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
shed = 2.0*(0.0082)*1.99 # similar to division but just a little bit different according to the plos epithelial strat paper

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


# Integration
   
print('integrating explicit')
x_path_poisson = odeintw(G, x_poisson, t)



print("Calculating")
 
for ti in range(0, time_total):
    cumu_extinct_poisson[ti] = np.sum(x_path_poisson[ti][0][:])
    


print('saving')
sys.stdout.flush()

with open('/gpfs1/home/m/c/mcboudre/scratch/hpv_modeling/cumu_extinct_2layer_main_poisson_11_27.npy', 'wb') as handle:
    np.save(handle, cumu_extinct_poisson)



print('done')
# sys.stdout.flush()