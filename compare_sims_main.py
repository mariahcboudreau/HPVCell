import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math


######
##
## SIMULATIONS
##
######

def cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t):
    # Set event times
    if skin[cell]['type'] == 'basal':
        tau_bbb = np.random.exponential(scale= 1/beta)
        tau_bbp = np.random.exponential(scale= 1/gamma)
        if tau_bbb < tau_bbp:
            tau_bpp = np.random.exponential(scale= 1/delta)
            if tau_bbb < tau_bpp:
                heapq.heappush(event_queue, (t + tau_bbb, cell, 'division - bbb'))
            else:
                heapq.heappush(event_queue, (t + tau_bpp, cell, 'division - bpp'))
        else:
            tau_bpp = np.random.exponential(scale=1 / delta)
            if tau_bbp < tau_bpp:
                heapq.heappush(event_queue, (t + tau_bbp, cell, 'division - bbp'))
            else:
                heapq.heappush(event_queue, (t + tau_bpp, cell, 'division - bpp'))


    elif skin[cell]['type'] == 'parabasal':
        tau_ppp = np.random.exponential(scale= 1/rho)
        tau_shed = np.random.exponential(scale= 1/theta)

        if tau_ppp < tau_shed:
            heapq.heappush(event_queue, (t + tau_ppp, cell, 'division - ppp'))
        else:
            heapq.heappush(event_queue, (t + tau_shed, cell, 'shed'))

    return event_queue

# Number of simulations to run
num_sims = 1000
all_history = []
all_basal_history = []
all_parabasal_history = []
all_times = []
all_bbb_count_history = []
all_bpp_count_history = []
all_bbp_count_history = []
all_ppp_count_history = []
all_shed_count_history = []

# Parameters of the model
symm_div = 0.04 
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
shed = 2.0*(0.0082)*1.99 # similar to division but just a little bit different according to the plos epithelial strat paper

rho = symm_div
gamma = asymm_div
delta = symm_div
beta = symm_div
theta = shed
basal_prop =  0.2 # TODO: Confirm proportion of basal cells
parabasal_prop = 1 - basal_prop
skin_size = 3*10**(6)

# Simulation
for s in range(num_sims):

    # Event queue for continuous time events
    event_queue = []

    # Queue for infected cells in the system
    skin = []

    # Initial conditions
    t = 0
    basals = 1
    parabasals = 0
    active_cells = parabasals + basals
    skin.append({'type': 'basal'})
    event_queue = cell_division(event_queue, skin, 0, beta, gamma, rho, delta, theta, t)

    # Counter

    history = []
    basal_history = []
    parabasal_history = []
    times = []
    bbb_count_history = []
    bpp_count_history = []
    bbp_count_history = []
    ppp_count_history = []
    shed_count_history = []
    history.append(active_cells / skin_size)
    basal_history.append(basals / skin_size)
    parabasal_history.append(parabasals / skin_size)
    bbb_count = 0
    bpp_count = 0
    bbp_count = 0
    ppp_count = 0
    shed_count = 0
    times.append(t)
    bbb_count_history.append(bbb_count)
    bpp_count_history.append(bpp_count)
    bbp_count_history.append(bbp_count)
    ppp_count_history.append(ppp_count)
    shed_count_history.append(shed_count)
   

    tmax = 1000
    while t < tmax and len(event_queue) > 0:
        (time, cell, event) = heapq.heappop(event_queue)
        t = time
        if active_cells < skin_size:
            if event == 'division - bbb' and (basals/skin_size) < basal_prop :
                skin.append({'type': 'basal'})
                event_queue = cell_division(event_queue, skin, len(skin)-1, beta, gamma,  rho, delta, theta, t)
                event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
                basals += 1
                bbb_count += 1
            elif event == 'division - bpp':
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                skin[cell]['type'] = 'dead'
                parabasals += 2
                basals -= 1
                bpp_count += 1
            elif event == "division - bbp":
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
                parabasals += 1
                bbp_count += 1
            elif event == "division - ppp":
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                parabasals += 1
                ppp_count += 1
            elif event == 'shed':
                skin[cell]['type'] = 'dead'
                parabasals -= 1
                shed_count += 1
        else:
            if event == "shed":
                skin[cell]['type'] = 'dead'
                parabasals -= 1
                shed_count += 1
        active_cells = parabasals + basals
        history.append(active_cells/skin_size)
        basal_history.append(basals/skin_size)
        parabasal_history.append(parabasals/skin_size)
        times.append(t)
        bbb_count_history.append(bbb_count)
        bpp_count_history.append(bpp_count)
        bbp_count_history.append(bbp_count)
        ppp_count_history.append(ppp_count)
        shed_count_history.append(shed_count)
    all_times.append(times)
    all_history.append(history)
    all_basal_history.append(basal_history)
    all_parabasal_history.append(parabasal_history)
    # all_bbb_count_history.append(bbb_count_history)
    # all_bpp_count_history.append(bpp_count_history)
    # all_bbp_count_history.append(bbp_count_history)
    # all_ppp_count_history.append(ppp_count_history)
    # all_shed_count_history.append(shed_count_history)



## Finding the average specifc time points. 
snapshots = [600, 800, 900]
avg_basals_arr = np.ndarray(shape = (len(snapshots),num_sims))
avg_parabasals_arr = np.ndarray(shape = (len(snapshots), num_sims))
avg_basals = np.zeros((len(snapshots)))
avg_parabasals = np.zeros((len(snapshots)))


for i in range(len(snapshots)):
    for sims in range(num_sims):
        temp_xb = np.linspace(0, int(all_times[sims][-1]), num = int(all_times[sims][-1])+1)
        temp_yb = np.interp(temp_xb, all_times[sims], all_basal_history[sims])
        if len(temp_xb) > snapshots[i]:
            avg_basals_arr[i][sims] = temp_yb[snapshots[i]]
        else: 
            avg_basals_arr[i][sims] = np.nan
        temp_xp = np.linspace(0, int(all_times[sims][-1]), num=int(all_times[sims][-1])+1)
        temp_yp = np.interp(temp_xp, all_times[sims], all_parabasal_history[sims])
        if len(temp_xp) > snapshots[i]:
            avg_parabasals_arr[i][sims] = temp_yp[snapshots[i]]
        else: 
            avg_parabasals_arr[i][sims] = np.nan
        
    
    avg_basals[i] = np.nanmean(avg_basals_arr[i][:])
    avg_parabasals[i] = np.nanmean(avg_parabasals_arr[i][:])



    


# fig, axs = plt.subplots(7, figsize=(6, 9), sharex= True)

# for i in range(num_sims-1):
#     axs[0].plot(all_times[i], all_bbb_count_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     axs[1].plot(all_times[i], all_bpp_count_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     axs[2].plot(all_times[i], all_bbp_count_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     axs[3].plot(all_times[i], all_basal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     axs[4].plot(all_times[i], all_ppp_count_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     axs[5].plot(all_times[i], all_parabasal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1,  rasterized=True)
#     axs[6].plot(all_times[i], all_shed_count_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)

# axs[0].set_ylabel('BBB Division')
# axs[1].set_ylabel('BPP Division')
# axs[2].set_ylabel('BBP Division')
# axs[3].set_ylabel('Basals')
# axs[4].set_ylabel('PPP Division')
# axs[5].set_ylabel('Parabasals')
# axs[6].set_ylabel('Shed')
# axs[6].set_xlabel('Time')



# plt.show()



#####
##
## MATHEMATICAL
##
#####

from numba import jit
from odeintw import odeintw

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

# Time of observations
t_length = 1000
t_steps = 1000
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions
nb_of_states = 100
nb_of_states_b = 20
nb_of_states_p = 80
x_0 = np.zeros((nb_of_states_b, nb_of_states_p))
x_0[1,0] = 1
m_0 = np.zeros(5)
m_0[0] = 1
m_0[2] = 1

# Parameters of the model
symm_div = 0.04
asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
shed = 2.0*(0.0082)*1.99 # similar to division but just a little bit different according to the plos epithelial strat paper

rho = symm_div
gamma = asymm_div
delta = symm_div
beta = symm_div
theta = shed



# Integration
# G = lambda x, t: J(x, t, beta, gamma, delta, rho, theta)
# x_path = odeintw(G, x_0, t_vec)
M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
m_path = odeintw(M, m_0, t_vec)



######
##
## PLOTTING
##
######

# plt.vlines(m_path[1][1], 0, 0.3, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[1][1]+stddev_parabasal, 0, 0.175, colors='gray', linestyles='--', label='Standard deviation')
# plt.vlines(m_path[1][1]-stddev_parabasal, 0, 0.325, colors='gray', linestyles='--', label='Standard deviation')


fig, axs = plt.subplots(3, sharex= True)

for i in range(num_sims-1):
    axs[0].plot(all_times[i], all_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
    axs[1].plot(all_times[i], all_basal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
    axs[2].plot(all_times[i], all_parabasal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
for t, j in zip(snapshots, range(len(snapshots))):
    std = np.sqrt(m_path[t][3] - m_path[t][1] ** 2)

    axs[2].plot(t, m_path[t][1]/skin_size, marker = "x", color = 'red')
    axs[2].plot(t, m_path[t][1]/skin_size + 2 * std, marker="_", color="red")
    axs[2].plot(t, avg_parabasals[j], marker = 'x', color = "yellow")
    # axs[2].plot(t, m_path[t][1]/skin_size - 2 * std, marker="_", color="red")

    std = np.sqrt(m_path[t][2]-m_path[t][0]**2)

    axs[1].plot(t, m_path[t][0]/skin_size, marker = "*", color = 'red')
    axs[1].plot(t, m_path[t][0]/skin_size + 2*std, marker = "_", color = "red")
    axs[1].plot(t, avg_basals[j], marker = 'x', color = "yellow")
    # axs[1].plot(t, m_path[t][0]/skin_size - 2 * std, marker="_", color="red")

axs[0].set_ylabel('Infected Cells')
axs[1].set_ylabel('Infected Basal Cells')
axs[2].set_ylabel('Infected Parabasal Cells')
axs[2].set_xlabel('Time')
axs[1].set_ylim([0,0.01])
axs[2].set_ylim([0,0.2])

plt.legend()
plt.show()
# plt.close()
# plt.savefig("sims_compare_mom.pdf", format = "pdf")
# print("done")

### 3D PLOTTING OF THE SIMS, ONLY A FEW SETS OF THEM
import matplotlib.ticker as ticker

# ax = plt.figure().add_subplot(projection='3d')
# majors = []


# for sim in range(0,num_sims,100):
#     ax.plot(all_times[sim], all_basal_history[sim], sim,  marker = "o", ls = '--', label = 'Simulation', rasterized=True)
#     majors.append(sim)
# for t, s, in zip(snapshots, range(len(snapshots))):
#     ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
#     ax.set_xlabel('Time')
#     ax.set_ylabel('Proportion of infected basal cells')
#     ax.set_zlabel('Simulation')
#     ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
# ax.stem(, y_vals, z_vals, linefmt='black', markerfmt = 'black', label='First moment')

# plt.show()

ax = plt.figure().add_subplot(projection='3d')
majors = []


for sim in range(0,num_sims,100):
    ax.plot(all_times[sim], all_parabasal_history[sim], sim,  marker = "o", ls = '--',  label = 'Simulation', rasterized=True)
    majors.append(sim)
ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
ax.set_xlabel('Time')
ax.set_ylabel('Proportion of infected parabasal cells')
ax.set_zlabel('Simulation')
ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
plt.show()