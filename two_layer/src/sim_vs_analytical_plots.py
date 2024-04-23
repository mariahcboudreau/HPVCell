import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math
import csv
import scipy



###################### PLOTTING ##############################

# Input values
num_sims = 50000
tmax = 500
dt=100


###################### OPENING SIMULATION FILES ##############################


from datetime import datetime
date = datetime.today().strftime('%m-%d-%Y')

date = "04-17-2024"

# with open('extinction_times_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     extinction_times = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/otherbasal_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    two_basal_history_other = np.load(f) 

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/basal_history_nonextinct_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    two_basal_history_nonextinct = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/dead_history_nonextinct_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    two_dead_history_nonextinct = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/otherparabasal_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    two_parabasal_history_other = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/otherdead_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    two_dead_count = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/otherextinct_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    two_extincts = np.load(f)

date = "03-15-2024"


# with open('extinction_times_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     extinction_times = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/otherbasal_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    four_basal_history_other = np.load(f) 

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/otherparabasal_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    four_parabasal_history_other = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/intermed_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    four_intermed_history = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/super_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    four_super_history = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/otherdead_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    four_dead_count = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/otherextinct_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    four_extincts = np.load(f)



###################### OPENING MOM FILES ##############################



# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'rb') as handle:    
#      extinct_mom_b_geometric = np.load(handle)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'rb') as f:
    four_extinction_mom_b_withdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/shed_first_moments_delta_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    four_shed_first_momentwithdead = np.load(f)

with open('four_layer/data/basal_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    four_basal_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/para_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
   four_para_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/intermed_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    four_intermed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/super_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    four_super_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/para_second_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    four_para_second_momentwithdead = np.load(f)



    
# date = "03-04-2024"

# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/extinction_mom_b_2layer_geometric_'+date+'_time500.npy', 'rb') as handle:    
#      extinct_mom_b_geometric = np.load(handle)


date = "04-10-2024"

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time500.npy', 'rb') as f:
    two_extinction_mom_b_withdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/shed_first_moments_delta_2layerwithdead_'+date+'_time500.npy', 'rb') as f:
    two_shed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/basal_first_moment_geom_2layerwithdead_'+date+'_time500.npy', 'rb') as f:
    two_basal_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/basal_second_moment_geom_2layerwithdead_'+date+'_time500.npy', 'rb') as f:
    two_basal_second_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/para_first_moment_geom_2layerwithdead_'+date+'_time500.npy', 'rb') as f:
    two_para_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/para_second_moment_geom_2layerwithdead_'+date+'_time500.npy', 'rb') as f:
    two_para_second_momentwithdead = np.load(f)


############################### COLORS ##################################
    

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


############################### PROBABILITY OF EXTINCTION ###########################################



plt.plot(np.arange(0,len(two_extincts)-1)/dt, two_extincts[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dashdot', linewidth = 3, label= "Simulated Three-agent - %d simulations" %(num_sims))
plt.plot(np.linspace(0, tmax, tmax), two_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments Three-agent')
plt.plot(np.arange(0,len(four_extincts)-1)/dt, four_extincts[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dashdot', linewidth = 3, label= "Simulated Five-agent - %d simulations" %(num_sims))
plt.plot(np.linspace(0, tmax, tmax), four_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments Five-agent')
plt.title("Probability of Extinction of Infected Basal Cells")
plt.xlabel('Time')
plt.ylabel('Probability')
plt.legend()
plt.show()
#plt.savefig('figures/mom-sims/prob_of_extinction_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()


############################### AVERAGE BASAL HISTORY - TWO ###########################################
date = datetime.today().strftime('%m-%d-%Y')

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.plot(np.arange(0,len(two_basal_history_other)-1)/dt, two_basal_history_other[1:]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax1.plot(np.linspace(0,tmax, num= tmax), two_basal_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
ax1.legend(loc= 'best')
ax1.set_ylabel("Basal")

############################### AVERAGE PARABASAL HISTORY ###########################################


ax2.plot(np.arange(0,len(two_parabasal_history_other))/dt, two_parabasal_history_other/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax2.plot(np.linspace(0,tmax, tmax), two_para_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
ax2.set_ylabel("Parabasal")


############################### AVERAGE DEAD HISTORY ###########################################


ax3.plot(np.arange(0,len(two_dead_count))/dt, two_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax3.plot(np.linspace(0,tmax, tmax), two_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
ax3.set_ylabel("Dead")

ax1.set_title("Averages of Agents")
#fig.savefig('figures/mom-sims/avg_basal_and_parabasal_dead_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
# fig.close()
############################### VARIANCE PARABASAL HISTORY ###########################################


# variance_mom = para_second_momentwithdead - para_first_momentwithdead**2
# vals = parabasal_variable
# vars = np.var(vals, axis = 0)

# plt.plot(np.arange(0,len(vars))/dt, vars, label = "Mid Sim variance")
# plt.plot(np.linspace(0,tmax, tmax), variance_mom, 'o', label = "MoM second moment estimations")
# plt.legend(loc= 'best')
# plt.title("Variance Parabasal Cells over Time - 2 Layer with Dead cells")
# plt.xlabel("Time")
# plt.ylabel("Variance Parabasal Cells")
# plt.show()
# plt.close()


############################### AVERAGE VIRIONS ################################################

plt.plot(np.arange(0,len(two_dead_count))/dt, 1000*two_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dashdot', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0, tmax, tmax), 1000*two_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
plt.title('Average Virions over Time')
plt.xlabel('Time')
plt.ylabel('Average Virons')
plt.legend()
plt.show()
#plt.savefig('figures/mom-sims/avg_virons_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()

############################### AVERAGE BASAL FOR PERSISTENT INFECTIONS ################################################

# Sims

two_basal_history_nonextinct[two_basal_history_nonextinct == 0] = np.nan
total_non_extincts = np.nanmean(two_basal_history_nonextinct, axis = 0)


# Analytical
mom_nonextincts = np.zeros((tmax))
for ti in range(0, tmax):
    p = ((2*two_basal_first_momentwithdead[ti])/(two_basal_first_momentwithdead[ti] + two_basal_second_momentwithdead[ti]))
    mom_nonextincts[ti] = 1/p


plt.plot(np.arange(0,len(total_non_extincts)-1)/dt, total_non_extincts[1:], c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0,tmax-1, num = tmax-1), mom_nonextincts[1:], c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
#plt.fill_between(np.arange(0,len(total_non_extincts_100)-1)/10, total_non_extincts_100[1:]+ var_total_non_extincts_100[1:], total_non_extincts_100[1:], alpha = 0.5)
plt.legend(loc= 'best')
plt.xlabel('Time')
plt.ylabel('Basals excluding extinction')
plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/figures/mom-sims/avg_nonextinct_basals_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()