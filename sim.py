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
num_sims = 1000
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
symm_div = 0.16    # Symmetric division rate - Clayton
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
basal_prop =  0.2 # TODO: Confirm proportion of basal cells 10^3 from plos epithelial
parabasal_prop = 1 - basal_prop
skin_size = 100000 #
# shed_amount = 100

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
  


    history.append(active_cells / skin_size)
    basal_history.append(basals / skin_size)
    parabasal_history.append(parabasals / skin_size)
    times.append(t)
    
  

    tmax = 500
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
            elif event == "shed":
                skin[cell]['type'] = 'dead'
                parabasals -= 1
                dead += 1
                all_shed_times.append(t)
                sheds.append(t)
                # dead_history.append(dead)
        else:
            if event == "shed":
                skin[cell]['type'] = 'dead'
                parabasals -= 1
                dead += 1
                all_shed_times.append(t)
                sheds.append(t)
                # dead_history.append(dead)
        active_cells = parabasals + basals
        history.append(active_cells/skin_size)
        basal_history.append(basals/skin_size)
        parabasal_history.append(parabasals/skin_size)
        times.append(t)
    all_times.append(times)
    all_history.append(history)
    all_basal_history.append(basal_history)
    all_parabasal_history.append(parabasal_history)
    shed_times_history.append(sheds)
    # all_dead.append(dead_history)
    # all_shed_times.append(shed_times)
    try: 
        ext_time = basal_history.index(0)
        extinction_times.append(times[basal_history.index(0)])
    except ValueError:
        ext_time = tmax
        extinction_times.append(tmax)



### Extinction probability computation

# extinction_count = np.zeros((tmax))

# for i in range(num_sims-1):
#     try:
#         index  = all_basal_history[i].index(0)
#     except ValueError:
#         index = len(all_times[i]) - 1
#     if np.ceil(all_times[i][index]) >= 200:
#         time_point = 199
#     else:
#         time_point = np.ceil(all_times[i][index])
#     if(time_point == 0):
#         print('stop here')
#     extinction_count[int(time_point)] += 1
        

with open('extinct_mom_b_2layer_main_geometric_1_8_500.npy', 'rb') as handle:    
     extinct_mom_b_geometric = np.load(handle)

# with open('cumu_extinct_2layer_main_delta_1_8.npy', 'rb') as f:
#     cumu_extinct_delta = np.load(f)

with open('shed_first_moments_delta_1-8_500.npy', 'rb') as f:
    shed_first_moment_delta = np.load(f)

# with open('shed_first_moments_poisson_12-11.npy', 'rb') as f:
#     shed_first_moment_poisson = np.load(f)

ext_count = np.zeros((tmax))
for i in range(0,tmax): 
    ext_count[i] = sum(j < i for j in extinction_times)

prob_extinction = ext_count/num_sims

dead_count = np.zeros((tmax))
ind_dead_count = np.zeros((num_sims, tmax))
for i in range(0,tmax):
    dead_count[i] = sum(j < i for j in all_shed_times)
    for s in range(num_sims):
        ind_dead_count[s,i] = sum(j < i for j in shed_times_history[s])

vars = np.zeros(tmax)
for i in range(0,tmax):
    vars[i] = np.var(ind_dead_count[:,i])

vars = 1000**2 * vars

rate_shed = dead_count/num_sims






plt.plot(np.linspace(0,tmax, tmax), prob_extinction, label= "Simulated probability of extinction - %d simulations" %(num_sims))
# plt.plot(np.linspace(0,200,200), cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta intital conditions)')
plt.plot(np.linspace(0, tmax, tmax), extinct_mom_b_geometric, label = 'Cumulative probability of extinction - MOM basals (Geometric approx & Delta Initial conditions)')
plt.title("Probabilty of Extinction - Simulations")
plt.xlabel('Time')
plt.ylabel('Probability')
plt.legend()
# plt.show()
plt.savefig('probextinct-sims1000-time500_1_8.pdf', format = 'pdf')
plt.close()




# discrete_death_array = np.zeros((tmax, num_sims))
# for t in range(1, tmax):
#     for i in range(num_sims-1):
#         index = [ti for ti in range(len(all_shed_times[i])) if all_shed_times[i][ti] < t]
#          # index of the last dead event time under discrete time point
#         if len(index) == 0:
#             discrete_death_array[t][i] = 0 
#         else:
#             discrete_death_array[t][i] = all_dead[i][index[-1]] 
        

# avg_dead_array = np.mean(discrete_death_array, axis = 1)
avg_virions_array = np.zeros((tmax))

avg_virions_array = 1000*rate_shed



plt.plot(np.linspace(0,tmax,tmax), avg_virions_array, label = "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0,tmax,tmax), shed_first_moment_delta, label = "MoM derivation - Delta initial conditions")
# for time in range(0,499, 50):
#     plt.errorbar(time, avg_virions_array[time], yerr = vars[time], ecolor= 'blue')
#eplt.plot(np.linspace(0,tmax,tmax), shed_first_moment_poisson, label = "MOM derivation - Poisson initial conditions")
plt.title("Average virions over time")
plt.xlabel('Time')
# plt.ylim(0, 25000)
plt.legend()
#plt.show()
# plt.savefig('avg_viral_load-sims1000-time500_1_8.pdf', format = 'pdf')


# Plotting flat simulations



fig, axs = plt.subplots(3)
#
for i in range(num_sims-1):

    # axs.plot(all_times[i], all_sheds[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, rasterized=True)
    # axs.set_ylabel('Virons shed')
    # axs.set_xlabel('Time')
    axs[0].plot(all_times[i], all_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation', rasterized=True)

    axs[1].plot(all_times[i], all_basal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation', rasterized=True)

    axs[2].plot(all_times[i], all_parabasal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation', rasterized=True)


axs[0].set_ylabel('Infected Cells')
axs[1].set_ylabel('Infected Basal Cells')
axs[2].set_ylabel('Infected Parabasal Cells')
axs[2].set_xlabel('Time')
plt.show()
# plt.savefig("~/Documents/MOCS2/testCode/HPVCellSim/shedding_simulation.pdf", format = 'pdf')


### 3D PLOTTING OF THE SIMS, ONLY A FEW SETS OF THEM
import matplotlib.ticker as ticker

ax = plt.figure().add_subplot(projection='3d')
majors = []


for sim in range(0,num_sims,100):
    ax.plot(all_times[sim], all_basal_history[sim], sim,  marker = "o", ls = '--', label = 'Simulation', rasterized=True)
    majors.append(sim)
   
ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
ax.set_xlabel('Time')
ax.set_ylabel('Proportion of infected basal cells')
ax.set_zlabel('Simulation')
ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
plt.show()

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