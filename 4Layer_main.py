# HPVCellSim Master Equations - 4 Layer

import numpy as np
import matplotlib.pyplot as plt
from numba import jit



# WE will use the odeint routine
from scipy.integrate import odeint
# With a wrapper to facilitate 2d arrays
from odeintw import odeintw

# Master Equations


@jit(nopython=True)
def J(c, t, beta, gamma, delta, rho, alpha, sigma, theta):
    """
        Time derivative of the occupation numbers.

            * c is the state distribution (array like)
            * t is time (scalar)
            * beta is the basal to basal/basal division rate
            * gamma is the basal to basal/parabasal division rate
            * delta is the basal to parabasal/parabasal division rate
            * rho is the parabasal to parabasal/parabasal divivsion rate
            * alpha is the parabasal transitioning to intermediate layer rate
            * sigma is the intermediate cells transitioning to the superficial layer rate
            * theta is the superficial layer shedding/death rate
            * NB: We use logistic growth for preys to limit the # of states
            * K will be the carrying capacity
        """
    dx = 0 * c
    


   
    # dx[b, p] is the entire master equation
    for b, p, i, s in np.ndindex(c.shape):
        if b < (c.shape[0] - 1): # basal to basal/basal output
            dx[b, p, i, s] -= beta * b * c[b, p, i, s]
        if p < c.shape[1] - 1:  # basal to parabasal/basal output
            dx[b, p, i, s] -= gamma * b * c[b, p, i, s]
        if p < c.shape[1] - 1:  # parabasal to parabasal/parabasal output
            dx[b, p, i, s] -= rho * p * c[b, p, i, s]
        if p < c.shape[1] - 2:  # basal to parabasal/parabasal output
            dx[b, p, i, s] -= delta * b * c[b, p, i, s]
        if p < c.shape[1] - 1: # parabasal to intermediate layer output
            dx[b, p, i, s] -= alpha * p * c[b, p, i, s]
        if i < c.shape[2] - 1: # intermediate to superficial layer output
            dx[b, p, i, s] -= sigma * i * c[b, p, i, s]
        if s < c.shape[3] - 1: # superficial shedding output
            dx[b, p, i, s] -= theta * s * c[b, p, i, s]
        if b > 1 :  # basal to basal/basal input
            dx[b, p, i, s] += beta * (b - 1) * c[b - 1, p, i, s]
        if b > 0 and p > 0:  # basal to parabasal/basal input
            dx[b, p, i, s] += gamma * b * c[b, p - 1, i, s]
        if p > 1:  # parabasal to parabasal/parabasal input
            dx[b, p, i, s] += rho * (p - 1) * c[b, p - 1, i, s]
        if b < (c.shape[0] - 1) and p > 1:  # basal to parabasal/parabasal input
            dx[b, p, i, s] += delta * (b + 1) * c[b + 1, p - 2, i, s]
        if p > 0: # parabasal to intermediate layer input
            dx[b, p, i, s] += alpha * (p + 1) * c[b, p + 1, i, s]
        if i > 0: # intermediate to superficial layer input
            dx[b, p, i, s] += sigma * (i + 1) * c[b, p, i + 1, s]
        if s < c.shape[1] - 1: # superficial shedding input
            dx[b, p, i, s] += theta * (s + 1) * c[b, p, i, s + 1]

    return dx




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
    dx[3] = -sigma * m[3]

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
    dx[10] = (beta - delta) * m[10]

    # Covariance for pi
    dx[11] = (rho - alpha - sigma) * m[11] + (gamma + 2 * delta) * m[9]

    # Covariance for ps
    dx[12] = (rho - theta - alpha) * m[12] + (gamma + 2 * delta) * m[10]

    # Covariance for is
    dx[13] = -1*(theta + sigma) * m[13]

    return dx

# Time of observations
t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions

nb_of_states_b = 10
nb_of_states_p = 10
nb_of_states_i = 10
nb_of_states_s = 10
x_0 = np.zeros((nb_of_states_b, nb_of_states_p, nb_of_states_i, nb_of_states_s))
x_0[1,0,0,0] = 1
m_0 = np.zeros(14)
m_0[0] = 1

# Parameters of the model
symm_div = 0.04
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
diff = 2.0*(0.0082)# similar to division but just a little bit different according to the plos epithelial strat paper

rho = symm_div
gamma = asymm_div
delta = symm_div
beta = symm_div
theta = 0.01 # plost epithelial chart
alpha = diff
sigma = 0.18 # plos epithelial chart 



# Integration
G = lambda x, t: J(x, t, beta, gamma, delta, rho, alpha, sigma, theta)
x_path = odeintw(G, x_0, t_vec)
M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, alpha, sigma, theta)
m_path = odeintw(M, m_0, t_vec)

print('done')
# extinct = np.zeros((t_length))
# cumu_extinct = np.zeros((t_length))
# for t in range(t_length - 1):
#     extinct[t] = x_path[t+1][0][0][0][0] - x_path[t][0][0][0][0]
#     cumu_extinct[t] = x_path[t][0][0][0][0]



# Axis 0 compresses the rows, Axis 1 compresses the columns
# AXIS 0 COMPRESSES THE FIRST ONE, AXIS 1 COMPRESSES THE SECOND


import matplotlib.ticker as ticker
ax = plt.figure().add_subplot(projection='3d')
majors = []
m_path_vals = []
y_vals = []
z_vals = []
stddev = np.sqrt(m_path[-1][4]-m_path[-1][0]**2)
for t in range(0,t_length-1,100):
    ax.plot(range(nb_of_states_b), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 1), axis = 1), t, marker="o", ls='--')
    majors.append(t)
    m_path_vals.append(m_path[t][0])
    y_vals.append(0)
    z_vals.append(t)
ax.stem(m_path_vals, y_vals, z_vals, linefmt='black', markerfmt = 'black', label='First moment')
ax.stem(m_path_vals+stddev, y_vals, z_vals, linefmt='gray', markerfmt = 'gray',  label='Standard deviation')
ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
ax.set_ylabel('Probabtility')
ax.set_xlabel('Number of infected parabasal cells')
ax.set_zlabel('Time')
ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
plt.show()
plt.close()

ax = plt.figure().add_subplot(projection='3d')
majors = []
m_path_vals = []
y_vals = []
z_vals = []
stddev = np.sqrt(m_path[-1][5]-m_path[-1][1]**2)
for t in range(0,t_length-1,100):
    ax.plot(range(nb_of_states_p), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 1), axis = 1), t, marker="o", ls='--')
    majors.append(t)
    m_path_vals.append(m_path[t][0])
    y_vals.append(0)
    z_vals.append(t)
ax.stem(m_path_vals, y_vals, z_vals, linefmt='black', markerfmt = 'black', label='First moment')
ax.stem(m_path_vals+stddev, y_vals, z_vals, linefmt='gray', markerfmt = 'gray',  label='Standard deviation')
ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
ax.set_ylabel('Probabtility')
ax.set_xlabel('Number of infected parabasal cells')
ax.set_zlabel('Time')
ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
plt.show()
plt.close()

# # # # Plot
# for t in np.arange(1, t_length-1, 100):
#     plt.plot(range(nb_of_states_b), np.sum(np.sum(np.sum(x_path[t], axis = 1), axis = 1), axis = 1), marker="o", ls='--',
#              label=fr"$t = {t_length * t / t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected basal cells')
# plt.show()
# #plt.savefig("explicit_basal.pdf", format = "pdf")
# plt.close()
# # # #
# # # # # Plot
# for t in np.arange(1,t_length-1,100):
#     plt.plot(range(nb_of_states_p), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 1), axis = 1), marker="o", ls='--', label=fr"$t = {t_length*t/t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected parabasal cells')
# plt.show()
# #plt.savefig("explicit_parabasal.pdf", format = "pdf")
# plt.close()

# # # # # Plot
# for t in np.arange(1,t_length-1,100):
#     plt.plot(range(nb_of_states_i), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 0), axis = 1), marker="o", ls='--', label=fr"$t = {t_length*t/t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected intermediate cells')
# plt.show()
# #plt.savefig("explicit_parabasal.pdf", format = "pdf")
# plt.close()

# # # # # Plot
# for t in np.arange(1,t_length-1,100):
#     plt.plot(range(nb_of_states_s), np.sum(np.sum(np.sum(x_path[t], axis = 1), axis = 1), axis = 1), marker="o", ls='--', label=fr"$t = {t_length*t/t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected superficial cells')
# plt.show()
# #plt.savefig("explicit_parabasal.pdf", format = "pdf")
# plt.close()

# stddev = np.sqrt(m_path[-1][4]-m_path[-1][0]**2)
# stddev_parabasal = np.sqrt(m_path[-1][5]-m_path[-1][1]**2)


# plt.vlines(m_path[10][0], 0, 0.9, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[10][0]+stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# #plt.vlines(m_path[1][0]-stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# plt.plot(range(nb_of_states_b), np.sum(np.sum(np.sum(x_path[10], axis = 1), axis = 1), axis = 1), marker="o", lw=0, ms=15, label=fr"$t = {t_vec[10]}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected basal cells')
# plt.show()
# plt.close()

# plt.vlines(m_path[10][1], 0, 0.9, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[10][1]+stddev_parabasal, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# #plt.vlines(m_path[1][0]-stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# plt.plot(range(nb_of_states_p), np.sum(np.sum(np.sum(x_path[10], axis = 0), axis = 1), axis = 1), marker="o", lw=0, ms=15, label=fr"$t = {t_vec[10]}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected parabasal cells')
# plt.show()
# plt.close()

# plt.plot(t_vec, extinct, label = 'Probability of extinction')
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()

# plt.plot(t_vec, cumu_extinct, label = 'Cumulative probability of extinction')
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()

