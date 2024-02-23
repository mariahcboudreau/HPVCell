import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math
import csv


####
#
#   Simulation Analyses
#
####


num_sims = 5000
tmax = 500
dt=100


######
# Opening simulated predictions
######

with open('extinction_times_02-22_time500_sims5000.npy', 'rb') as f:
    extinction_times = np.load(f)


# with open('all_times_02-22_time500_sims5000.csv', 'r') as read_obj:
#     csv_reader = csv.reader(read_obj)
#     all_times = list(csv_reader)
#     all_times = [list(map(float, sublist)) for sublist in all_times]

# with open('all_history_02-22_time500_sims5000.csv', 'r') as read_obj:
#     csv_reader = csv.reader(read_obj)
#     all_history = list(csv_reader)
#     all_history = [list(map(float, sublist)) for sublist in all_history]

# with open('all_basal_history_02-22_time500_sims5000.csv', 'r') as read_obj:
#     csv_reader = csv.reader(read_obj)
#     all_basal_history = list(csv_reader)
#     all_basal_history = [list(map(float, sublist)) for sublist in all_basal_history]

# with open('all_parabasal_history_02-22_time500_sims5000.csv', 'r') as read_obj:
#     csv_reader = csv.reader(read_obj)
#     all_parabasal_history = list(csv_reader)
#     all_parabasal_history = [list(map(float, sublist)) for sublist in all_parabasal_history]

# with open('shed_times_history_02-22_time500_sims5000.csv', 'r') as read_obj:
#     csv_reader = csv.reader(read_obj)
#     shed_times_history = list(csv_reader)
#     shed_times_history = [list(map(float, sublist)) for sublist in shed_times_history]

# with open('all_shed_times_02-22_time500_sims5000.npy', 'rb') as read_obj:
#     all_shed_times = np.load(read_obj)

with open('otherbasal_history_02-23_time500_sims5000.npy', 'rb') as f:
    basal_history_other = np.load(f) 

with open('otherparabasal_history_02-23_time500_sims5000.npy', 'rb') as f:
    parabasal_history_other = np.load(f)

with open('otherdead_history_02-23_time500_sims5000.npy', 'rb') as f:
    dead_count = np.load(f)


######
# Opening Analytical predictions
######
    
with open('extinction_mom_b_2layer_geometric_02-23_time500.npy', 'rb') as handle:    
     extinct_mom_b_geometric = np.load(handle)

# with open('cumu_extinct_2layer_main_delta_1_8.npy', 'rb') as f:
#     cumu_extinct_delta = np.load(f)

with open('shed_first_moments_delta_02-23_time500.npy', 'rb') as f:
    shed_first_moment_delta = np.load(f)

# with open('para_first_moments_geom_1-15_500.npy', 'rb') as f:
#     para_first_moment = np.load(f)

with open('basal_first_moment_geom_1-29_500.npy', 'rb') as f:
    basal_first_moment= np.load(f)

with open('para_first_moment_geom_1-29_500.npy', 'rb') as f:
    para_first_moment = np.load(f)

######
# Calculating probability of extinction
######


ext_count = np.zeros((tmax))
for i in range(0,tmax): 
    ext_count[i] = sum(j < i for j in extinction_times)

prob_extinction = ext_count/num_sims

######
# Calculating the rate of shed
######

rate_shed = dead_count/num_sims

avg_virions_array = 1000*rate_shed



plt.plot(np.arange(0,len(avg_virions_array))/dt, avg_virions_array, marker = '+', linewidth = 2, label = "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0,tmax,tmax), shed_first_moment_delta, label = "MoM derivation - Delta initial conditions")
# plt.plot(np.linspace(0,tmax,tmax), avg_virions_array + margin_error, c = 'grey')
# plt.plot(np.linspace(0,tmax,tmax), avg_virions_array - margin_error, c = 'grey', label = "95% CI")
# for time in range(0,499, 50):
#     plt.errorbar(time, avg_virions_array[time], yerr = vars[time], ecolor= 'blue')
plt.title("Average virions over time")
plt.xlabel('Time')
# plt.ylim(0, 25000)
plt.legend()
plt.show()
# plt.savefig('avg_viral_load-sims1000-time500_1_8.pdf', format = 'pdf')


# dead_count = np.zeros((tmax))
# ind_dead_count = np.zeros((num_sims, tmax))
# for i in range(0,tmax):
#     dead_count[i] = sum(j < i for j in all_shed_times)
#     for s in range(num_sims):
#         ind_dead_count[s,i] = sum(j < i for j in shed_times_history[s])

# vars = np.zeros(tmax)
# for i in range(0,tmax):
#     vars[i] = np.var(ind_dead_count[:,i])

# vars = 1000**2 * vars



# ######
# # Calculating margin of error or the standard deviation
# ###### 

# margin_error = np.zeros(tmax)
# for t in range(tmax):
#     margin_error[t] = 1.96 * ( np.sqrt(vars[t])/ np.sqrt(num_sims)) 


# ######
# # Interpolation of basal averages
# ###### 
    
# x_new = np.linspace(0,500, num= 500)
# y_new = np.zeros(([num_sims, 500]))
# for s in range(num_sims):
#     x_true = all_times[s]
#     y_true = all_basal_history[s]
#     y_new[s] = np.interp(x_new, x_true, y_true)


# y_new_avg = np.mean(y_new, axis = 0)


# plt.plot(x_new, y_new_avg, "-", label = "Linear interp from simulations")
plt.plot(np.arange(0,len(basal_history_other))/dt, basal_history_other/num_sims, label = "Mid Sim averages")
plt.plot(np.linspace(0,tmax, num= tmax), basal_first_moment, 'o', label = "MoM estimations")
plt.legend(loc= 'best')
plt.xlabel("Time")
plt.ylabel("Average Basal Cells")
plt.show()

# ######
# # Interpolation of parabasals averages
# ###### 

# x_new = np.linspace(0,500, num= 500)
# y_new = np.zeros(([num_sims, 500]))
# for s in range(num_sims):
#     x_true = all_times[s]
#     y_true = all_parabasal_history[s]
#     y_new[s] = np.interp(x_new, x_true, y_true)


# y_new_avg = np.mean(y_new, axis = 0)


# plt.plot(x_new, y_new_avg, "-", label = "Linear interp from simulations")
plt.plot(np.arange(0,len(parabasal_history_other))/dt, parabasal_history_other/num_sims, label = "Mid Sim averages")
plt.plot(np.linspace(0,tmax, tmax), para_first_moment, 'o', label = "MoM estimations")
plt.legend(loc= 'best')
plt.xlabel("Time")
plt.ylabel("Average Parabasal Cells")
plt.show()
plt.close()


# # 
# para_avg = np.zeros(tmax)
# for s in range(num_sims):
#     for i in range(len(all_times[s])-1):
#          for j in range(int(np.ceil(all_times[s])[i]), int(np.ceil(all_times[s])[i+1]-1)):
#                 if j >= 500:
#                     para_avg[-1] += all_parabasal_history[s][i]
#                 else:
#                     para_avg[j] += all_parabasal_history[s][i]

# para_avg = para_avg/num_sims

# plt.plot(np.linspace(0,tmax-1,tmax-1), para_first_moment[0:499], label = 'MoM para average')
# plt.plot(np.linspace(0,tmax-1,tmax-1), para_avg[0:499], label = 'Simulations para average')
# plt.legend()
# plt.show()
# # 
# b_avg = np.zeros(tmax)
# for s in range(num_sims):
#     for i in range(len(all_times[s])-1):
#          for j in range(int(np.ceil(all_times[s])[i]), int(np.ceil(all_times[s])[i+1]-1)):
#                 if j >= 500:
#                     b_avg[-1] += all_basal_history[s][i]
#                 else:
#                     b_avg[j] += all_basal_history[s][i]

# b_avg = b_avg/num_sims

# plt.plot(np.linspace(0,tmax-1,tmax-1), basal_first_moment[0:499], label = 'MoM b average')
# plt.plot(np.linspace(0,tmax-1,tmax-1), b_avg[0:499], label = 'Simulations b average')
# plt.legend()
# plt.show()


# plt.plot(np.linspace(0,tmax,tmax), shed_first_moment_delta, label = 'MoM Parabasal average times theta')
# plt.plot(np.linspace(0,tmax,tmax), rate_shed, label = 'Simulated rate of shed')
# plt.legend()
# plt.show()

plt.plot(np.linspace(0,tmax, tmax), prob_extinction, marker = "+", linewidth = 2, label= "Simulated probability of extinction - %d simulations" %(num_sims))
# plt.plot(np.linspace(0,200,200), cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta intital conditions)')
plt.plot(np.linspace(0, tmax, tmax), extinct_mom_b_geometric, linewidth = 2, label = 'Cumulative probability of extinction - MoM')
plt.title("Probabilty of Extinction of Infected Basal Cells")
plt.xlabel('Time')
plt.ylabel('Probability')
plt.legend()
plt.show()
#plt.savefig('probextinct-geom_sims100-time500_02_23.pdf', format = 'pdf')
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



# with open('simulated_probextinct_geom_numsim-million_infpop_date1-9_time500.npy', 'wb') as f:
#     np.save(f, prob_extinction)
# with open('simulated_avg-virions_geom_numsim-million_infpop_date1-9_time500.npy', 'wb') as f:
#     np.save(f, avg_virions_array)

# # Plotting flat simulations



# fig, axs = plt.subplots(3)
# #
# for i in range(num_sims-1):

#     # axs.plot(all_times[i], all_sheds[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, rasterized=True)
#     # axs.set_ylabel('Virons shed')
#     # axs.set_xlabel('Time')
#     axs[0].plot(all_times[i], all_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation', rasterized=True)

#     axs[1].plot(all_times[i], all_basal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation', rasterized=True)

#     axs[2].plot(all_times[i], all_parabasal_history[i], marker = "o", ls = '--', color = 'black', alpha = 0.25, label = 'Simulation', rasterized=True)


# axs[0].set_ylabel('Infected Cells')
# axs[1].set_ylabel('Infected Basal Cells')
# axs[2].set_ylabel('Infected Parabasal Cells')
# axs[2].set_xlabel('Time')
# plt.show()
# # plt.savefig("~/Documents/MOCS2/testCode/HPVCellSim/shedding_simulation.pdf", format = 'pdf')


# ### 3D PLOTTING OF THE SIMS, ONLY A FEW SETS OF THEM
# import matplotlib.ticker as ticker

# ax = plt.figure().add_subplot(projection='3d')
# majors = []


# for sim in range(0,num_sims,100):
#     ax.plot(all_times[sim], all_basal_history[sim], sim,  marker = "o", ls = '--', label = 'Simulation', rasterized=True)
#     majors.append(sim)
   
# ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
# ax.set_xlabel('Time')
# ax.set_ylabel('Proportion of infected basal cells')
# ax.set_zlabel('Simulation')
# ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
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