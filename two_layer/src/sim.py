# Author: Mariah Boudreau


################################# IMPORTS ###q################################
import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
from datetime import datetime
date = datetime.today().strftime('%m-%d-%Y')


####
#   cell_division function
#   Description: Defines the events occuring in the 2 Layer HPV cell division
#   
#   Inputs: event_queue - Events in order of which they will occur
#           skin        - List of infected cells with their type 
#           cell        - Index of cell having an event occur
#           beta        - Rate of bbb division
#           gamma       - Rate of bbp division
#           rho         - Rate of ppp division
#           delta       - Rate of bpp division
#           theta       - Rate of shed 
#           t           - Time
#
#   Outputs: Updated event queue with new event added
####
def cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t):

    # For basal cells, set event times
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
num_sims = 50000

# Next time differential
dt = 10

# Max time for the simulation to run
tmax = 750

# All history variables initialized
all_history = []
all_basal_history = []
all_parabasal_history = []
all_times = []
all_shed_times = []
extinction_times = []
shed_times_history = []


parabasal_history_other = np.zeros(tmax*dt)
parabasal_variable = np.zeros((num_sims,tmax*dt))
basal_history_other = np.zeros(tmax*dt)
basal_history_other[0] = num_sims
basal_history_nonextinct = np.zeros((num_sims,tmax*dt))
dead_history_other = np.zeros(tmax*dt)
dead_history_all = np.zeros((num_sims,tmax*dt))
extincts = np.zeros(tmax*dt)
extinct_sims = []

# Updated cell division parameters
R_b = 0.03          # Division rate of basal - Murall citation
R_p = 0.39          # Division rate of parabasal - Murall citation
symm_div = 0.08     # Symmetric division rate - Clayton
asymm_div = 0.84    # Asymmetric division rate - Clayton

beta = R_b * symm_div  + 0.001      # bbb with offset to avoid steady state
gamma = R_b * asymm_div             # bbp
delta = R_b * symm_div              # bpp
rho = R_p * symm_div                # ppp
theta = 0.67                        # Shed - Murall






################################# SIMULATION ###################################
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

    # Inital conditions for the infected cells
    skin.append({'type': 'basal'})
    event_queue = cell_division(event_queue, skin, 0, beta, gamma, rho, delta, theta, t)
    

    
    # Set time differential
    next_time = 1/dt

    # Set variable for extinction counter
    first = 0

    # Start timer
    while t < tmax:
        if len(event_queue) == 0:
            t = tmax
        else:
        # Pop the next event out
            (time, cell, event) = heapq.heappop(event_queue)
            t = time
         # Assign new values to the 
        while t > next_time and int(round(next_time*dt))<len(parabasal_history_other):
            # parabasal_history_other[int(round(next_time*dt))] += parabasals
            # # parabasal_variable[s, int(round(next_time*dt))] += parabasals
            # basal_history_other[int(round(next_time*dt))] += basals
            # basal_history_nonextinct[s, int(round(next_time*dt))] += basals
            # dead_history_other[int(round(next_time*dt))] += dead
            dead_history_all[s, int(round(next_time*dt))] += dead
            
            # Track the extinction events
            
            if basals == 0 and first == 0:
                extincts[int(round(next_time*dt)):-1] += 1
                first = 99
                extinct_sims.append(s)
            # Set next time value    
            next_time += 1/dt
        
        # bbb division
        if event == 'division - bbb':
            # Add new basal cell
            skin.append({'type': 'basal'})
            # Create new division events for the current cell and new cell
            event_queue = cell_division(event_queue, skin, len(skin)-1, beta, gamma, rho, delta, theta, t)
            event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
            # Add to the basal count
            basals += 1

        # bpp division
        elif event == 'division - bpp':
            # Add new parabasal cell
            skin.append({'type': 'parabasal'})
            # Create new division events for the new cell
            event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
            # Add new parabasal cell
            skin.append({'type': 'parabasal'})
            # Create new division events for the new cell
            event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
            # Reassign basal cell that died
            skin[cell]['type'] = 'dead basal'

            # Add to the parabasal count and take away from basal count
            parabasals += 2
            basals -= 1

        # bbp division
        elif event == "division - bbp":
            # Add new parabasal cell
            skin.append({'type': 'parabasal'})
            # Create new division events for the current cell and new cell
            event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
            event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
            # Add to the parabasal count
            parabasals += 1

        # ppp division
        elif event == "division - ppp":
            # Add new parabasal cell
            skin.append({'type': 'parabasal'})
            # Create new division events for the current cell and new cell
            event_queue = cell_division(event_queue, skin, len(skin) - 1, beta, gamma, rho, delta, theta, t)
            event_queue = cell_division(event_queue, skin, cell, beta, gamma, rho, delta, theta, t)
            # Add to the parabasal count
            parabasals += 1

        # shed 
        elif event == "shed":
            # Reassign parabasal cell that died
            skin[cell]['type'] = 'dead parabasal'
            # Change parabasal and dead counts
            parabasals -= 1
            dead += 1
            # all_shed_times.append(t)
            # sheds.append(t)
            # dead_history.append(dead)
        
  



############################## SAVING FILES #######################################


# with open('extinction_times_03-01_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
#     np.save(f, extinction_times)

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
    


# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/otherbasal_history_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
#     np.save(f, basal_history_other)

# with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherparabasal_history_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
#     np.save(f, parabasal_history_other)

basal_history_nonextinct[basal_history_nonextinct == 0] = np.nan
two_total_non_extincts = np.nanmean(basal_history_nonextinct, axis = 0)

two_dead_history_nonextinct = np.delete(dead_history_all, extinct_sims, axis = 0)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/basal_history_nonextinct_mean_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, two_total_non_extincts)

# with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherdead_history_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
#     np.save(f, dead_history_other)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/two_dead_history_nonextinct_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, two_dead_history_nonextinct)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/two_dead_history_nonextinct_mean_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, two_dead_history_nonextinct.mean(axis = 0))

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/two_dead_history_nonextinct_variance_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
    np.save(f, two_dead_history_nonextinct.var(axis = 0))



# with open('otherparabasal_variable_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
#     np.save(f, parabasal_variable)

# with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherextinct_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'wb') as f:
#     np.save(f, extincts)



############################## GET FUNCTIONS #######################################
    
def get_sim_basals():
    return basals

def get_sim_parabasals():
    return parabasals

def get_time():
    return t

def get_event():
    return event

def get_event_queue():
    return event_queue

def get_dead():
    return dead

def get_next_time():
    return next_time