import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math




def cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t):
    # Set event times
    if skin[cell]['type'] == 'basal':
        tau_bbb = np.random.exponential(scale= 1/beta)
        tau_bbp = np.random.exponential(scale= 1/gamma)
        if tau_bbb < tau_bbp:
            tau_bpp = np.random.exponential(scale=1 / delta)
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
num_sims = 100
all_history = []
all_basal_history = []
all_parabasal_history = []
all_times = []
all_sheds = []

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
skin_size = 1000
shed_amount = 100

# Simulation
for s in range(num_sims):

    # Event queue for continuous time events
    event_queue = []

    # Queue for infected cells in the system
    skin = []

    # Initial conditions
    shed = 0
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
    shed_history = []
    times = []
    history.append(active_cells / skin_size)
    basal_history.append(basals / skin_size)
    parabasal_history.append(parabasals / skin_size)
    times.append(t)
    shed_history.append(shed)

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
            elif event == 'division - bpp':
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                skin[cell]['type'] = 'dead'
                parabasals += 2
                basals -= 1
            elif event == "division - bbp":
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
                parabasals += 1
            elif event == "division - ppp":
                skin.append({'type': 'parabasal'})
                event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
                parabasals += 1
            elif event == 'shed':
                skin[cell]['type'] = 'dead'
                parabasals -= 1
                shed += shed_amount
        else:
            if event == "shed":
                skin[cell]['type'] = 'dead'
                parabasals -= 1
                shed += shed_amount
        active_cells = parabasals + basals
        history.append(active_cells/skin_size)
        basal_history.append(basals/skin_size)
        parabasal_history.append(parabasals/skin_size)
        shed_history.append(shed)
        times.append(t)
    all_times.append(times)
    all_history.append(history)
    all_basal_history.append(basal_history)
    all_parabasal_history.append(parabasal_history)
    all_sheds.append(shed_history)




# Plotting



fig, axs = plt.subplots(1)
#
for i in range(num_sims-1):

    axs.plot(all_times[i], all_sheds[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, rasterized=True)
    axs.set_ylabel('Virons shed')
    axs.set_xlabel('Time')
#     axs[0].plot(all_times[i], all_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation')
#
#     axs[1].plot(all_times[i], all_basal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation')
#
#     axs[2].plot(all_times[i], all_parabasal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation')
#
#
# axs[0].set_ylabel('Infected Cells')
# axs[1].set_ylabel('Infected Basal Cells')
# axs[2].set_ylabel('Infected Parabasal Cells')
# axs[2].set_xlabel('Time')
#plt.show()
plt.savefig("~/Documents/MOCS2/testCode/HPVCellSim/shedding_simulation.pdf", format = 'pdf')
