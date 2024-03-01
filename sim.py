import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math


####
#
#   
#
####
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
            tau_bpp = np.random.exponential(scale= 1/delta)
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
num_sims = 10000
all_history = []
all_basal_history = []
all_parabasal_history = []
all_times = []
all_shed_times = []
extinction_times = []
shed_times_history = []

# Updated parameters
R_b = 0.03          # Division rate of basal - Murall citation
R_p = 0.39          # Division rate of parabasal - Murall citation
symm_div = 0.08    # Symmetric division rate - Clayton
asymm_div = 0.84    # Asymmetric division rate - Clayton

beta = R_b * symm_div  + 0.001      # bbb 
gamma = R_b * asymm_div             # bbp
delta = R_b * symm_div              # bpp
rho = R_p * symm_div                # ppp
theta = 0.67                        # Shed - Murall


# # Parameters of the model - OLD

# # Probability of type of division
# symm_div = 0.04
# asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance

# # Rate of types of divisions


# # Rate of shedding
# shed = 2.0*(0.0082) # similar to division but just a little bit different according to the plos epithelial strat paper

# rho = symm_div
# gamma = asymm_div
# delta = symm_div
# beta = symm_div + 0.01
# theta = shed
basal_prop =  0.1 # TODO: Confirm proportion of basal cells 10^3 from plos epithelial
parabasal_prop = 1 - basal_prop
skin_size = 100000 #
# shed_amount = 100

dt = 100
tmax = 500
parabasal_history_other = np.zeros(tmax*dt)
parabasal_variable = np.zeros((num_sims,tmax*dt))
basal_history_other = np.zeros(tmax*dt)
basal_history_other[0] = 1
dead_history_other = np.zeros(tmax*dt)




# Simulation
for s in range(num_sims):

    # Event queue for continuous time events
    event_queue = []

    # Queue for infected cells in the system
    skin = []

    # Initial conditions
    dead = 0
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
    dead_history = []
    times = []
    sheds = []
  


    history.append(active_cells)
    basal_history.append(basals)
    parabasal_history.append(parabasals)
    times.append(t)
    

##
# Infinite population
##
   
    next_time = 1/dt
    while t < tmax and len(event_queue) > 0:
        (time, cell, event) = heapq.heappop(event_queue)
        t = time
        while time > next_time and int(round(next_time*dt))<len(parabasal_history_other):
            parabasal_history_other[int(round(next_time*dt))] += parabasals
            parabasal_variable[s, int(round(next_time*dt))] += parabasals
            basal_history_other[int(round(next_time*dt))] += basals
            dead_history_other[int(round(next_time*dt))] += dead
            next_time += 1/dt
        if event == 'division - bbb':
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
            event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
            parabasals += 1
        elif event == "shed":
            skin[cell]['type'] = 'dead'
            parabasals -= 1
            dead += 1
            all_shed_times.append(t)
            sheds.append(t)
            # dead_history.append(dead)
        
    #     active_cells = parabasals + basals
    #     history.append(active_cells)
    #     basal_history.append(basals)
    #     parabasal_history.append(parabasals)
    #     times.append(t)
    # all_times.append(times)
    # all_history.append(history)
    # all_basal_history.append(basal_history)
    # all_parabasal_history.append(parabasal_history)
    # shed_times_history.append(sheds)
        # all_dead.append(dead_history)
        # all_shed_times.append(shed_times)
    try: 
        ext_time = basal_history.index(0)
        extinction_times.append(times[basal_history.index(0)])
    except ValueError:
        ext_time = tmax
        extinction_times.append(tmax)






import csv

with open('extinction_times_03-01_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, extinction_times)

# with open("all_times_02-22_time500_sims1000.csv", "w") as f:
#     wr = csv.writer(f)
#     wr.writerows(all_times)   
 
# with open('all_history_02-22_time500_sims1000.csv', 'w') as f:
#     wr = csv.writer(f)
#     wr.writerows(all_history)

# with open('all_basal_history_02-22_time500_sims1000.csv', 'w') as f:
#     wr = csv.writer(f)
#     wr.writerows(all_basal_history)

# with open('all_parabasal_history_02-22_time500_sims1000.csv', 'w') as f:
#     wr = csv.writer(f)
#     wr.writerows(all_parabasal_history)

# with open('shed_times_history_02-22_time500_sims1000.csv', 'w') as f:
#     wr = csv.writer(f)
#     wr.writerows(shed_times_history)

# with open('all_shed_times_02-22_time500_sims1000.npy', 'wb') as f:
#     np.save(f,all_shed_times)

with open('otherbasal_history_03-01_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, basal_history_other)

with open('otherparabasal_history_03-01_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, parabasal_history_other)

with open('otherdead_history_03-01_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, dead_history_other)

with open('otherparabasal_variable_03-01_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, parabasal_variable)