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
snapshots = [100, 200, 300, 400, 475]
avg_basals_arr = np.ndarray(shape = (len(snapshots),num_sims))
avg_parabasals_arr = np.ndarray(shape = (len(snapshots), num_sims))
avg_basals = np.zeros((len(snapshots)))
avg_parabasals = np.zeros((len(snapshots)))
extinction = np.zeros((tmax))


for i in range(len(snapshots)):


    for sims in range(num_sims):
        temp_xb = np.linspace(0, int(all_times[sims][-1]), num = int(all_times[sims][-1])+1)
        temp_yb = np.interp(temp_xb, all_times[sims], all_basal_history[sims])
        
        extinction[np.where(temp_yb == 0)[0][0]] += 1
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
    extinction_prob = extinction/num_sims
    
    print(extinction_prob)


    


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



with open('m_path_2layer_main_1.npy', 'rb') as f:
    m_path = np.load(f)
with open('x_path_file_2layer_main.npy', 'rb') as f:
    x_path= np.load(f)


######
##
## PLOTTING
##
######

# plt.vlines(m_path[1][1], 0, 0.3, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[1][1]+stddev_parabasal, 0, 0.175, colors='gray', linestyles='--', label='Standard deviation')
# plt.vlines(m_path[1][1]-stddev_parabasal, 0, 0.325, colors='gray', linestyles='--', label='Standard deviation')


# fig, axs = plt.subplots(3, sharex= True)

# #for i in range(num_sims-1):
#     # axs[0].plot(all_times[i], all_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     # axs[1].plot(all_times[i], all_basal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
#     # axs[2].plot(all_times[i], all_parabasal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.1, rasterized=True)
# for t, j in zip(snapshots, range(len(snapshots))):
#     std = np.sqrt(m_path[t][3] - m_path[t][1] ** 2)

#     axs[2].plot(t, m_path[t][1]/skin_size, marker = "x", color = 'red')
#     #axs[2].plot(t, m_path[t][1]/skin_size + 2 * std, marker="_", color="red")
#     axs[2].plot(t, avg_parabasals[j], marker = 'x', color = "yellow")
#     # axs[2].plot(t, m_path[t][1]/skin_size - 2 * std, marker="_", color="red")

#     std = np.sqrt(m_path[t][2]-m_path[t][0]**2)

#     axs[1].plot(t, m_path[t][0]/skin_size, marker = "*", color = 'red')
#     markers1, = axs[1].plot(t, m_path[t][0]/skin_size + 2*std, marker = "_", color = "red")
#     markers2, = axs[1].plot(t, avg_basals[j], marker = 'x', color = "yellow")
#     if t == snapshots[0]:
#         markers1.set_label("Method of moments")
#         markers2.set_label("Average of simulations")
#     # axs[1].plot(t, m_path[t][0]/skin_size - 2 * std, marker="_", color="red")

# axs[0].set_ylabel('Infected Cells')
# axs[1].set_ylabel('Infected Basal Cells')
# axs[2].set_ylabel('Infected Parabasal Cells')
# axs[2].set_xlabel('Time')
# axs[1].set_ylim([0,0.005])
# axs[2].set_ylim([0,0.1])
# axs[1].legend()


# plt.show()

# #plt.savefig("sims_compare_mom.pdf", format = "pdf")
# plt.close()
# # print("done")

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

# ax = plt.figure().add_subplot(projection='3d')
# majors = []


# for sim in range(0,num_sims,100):
#     ax.plot(all_times[sim], all_parabasal_history[sim], sim,  marker = "o", ls = '--',  label = 'Simulation', rasterized=True)
#     majors.append(sim)
# ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
# ax.set_xlabel('Time')
# ax.set_ylabel('Proportion of infected parabasal cells')
# ax.set_zlabel('Simulation')
# ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
# plt.show()




#####
# Extinction probabilty
#####

cumu_extinct_delta = np.zeros((500))
extinct_mom_b_delta = np.zeros((500))

for time in range(500):
    cumu_extinct_delta[time] = np.sum(x_path[0][:])
    # MOM
    if m_path[0][2] != 0:
        extinct_mom_b_delta[time] = (1-((m_path[0][0]**2)/(m_path[0][2])))
    else: 
        extinct_mom_b_delta[time] = 0

t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)


plt.plot(t_vec, cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta Approx)')
#plt.plot(t_vec, extinct_mom_b_delta, label = 'Cumulative probability of extinction - MOM basals (Delta Approx)')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('2 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()


plt.plot(t_vec, extinct_mom_b_delta, label = 'Cumulative probability of extinction - MOM basals (Delta Approx)')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('2 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()