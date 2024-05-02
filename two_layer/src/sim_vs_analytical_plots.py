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
tmax = 750
dt=100


###################### OPENING SIMULATION FILES ##############################


from datetime import datetime

date = "04-22-2024"

# with open('extinction_times_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     extinction_times = np.load(f)

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherbasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_basal_history_other = np.load(f) 

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/basal_history_nonextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_basal_history_nonextinct = np.load(f)

# with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/dead_history_nonextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
#     two_dead_history_nonextinct = np.load(f)

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherparabasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_parabasal_history_other = np.load(f)

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherdead_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_dead_count = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/two_layer_data/otherextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_extincts = np.load(f)


##### Four Layer data 

# with open('extinction_times_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     extinction_times = np.load(f)

date = "04-23-2024"

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/four_layer_data/otherbasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_basal_history_other = np.load(f) 

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/four_layer_data/otherparabasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_parabasal_history_other = np.load(f)

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/four_layer_data/intermed_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_intermed_history = np.load(f)

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/four_layer_data/super_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_super_history = np.load(f)

with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/four_layer_data/otherdead_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_dead_count = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('/Users/mcboudre/OneDrive - University of Vermont/HPV-Data/four_layer_data/otherextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_extincts = np.load(f)



###################### OPENING MOM FILES ##############################

date = "04-25-2024"

# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'rb') as handle:    
#      extinct_mom_b_geometric = np.load(handle)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_extinction_mom_b_withdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/shed_first_moments_delta_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_shed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/basal_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_basal_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/para_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
   four_para_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/intermed_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_intermed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/super_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_super_first_momentwithdead = np.load(f)

# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/para_second_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
#     four_para_second_momentwithdead = np.load(f)




# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/extinction_mom_b_2layer_geometric_'+date+'_time500.npy', 'rb') as handle:    
#      extinct_mom_b_geometric = np.load(handle)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_extinction_mom_b_withdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/shed_first_moments_delta_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_shed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/basal_first_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_basal_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/basal_second_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_basal_second_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/para_first_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_para_first_momentwithdead = np.load(f)

# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/para_second_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
#     two_para_second_momentwithdead = np.load(f)


############################### COLORS ##################################
    

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


############################### PROBABILITY OF EXTINCTION - Both Layers ###########################################
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
plt.rcParams["font.family"] = "Times New Roman"

fig, (ax1, ax2) = plt.subplots(2, sharex=True)


ax1.plot(np.arange(0,len(two_extincts)-1)/dt, two_extincts[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated Three-agent")
ax1.plot(np.linspace(0, tmax, tmax), two_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments Three-agent')
ax2.plot(np.arange(0,len(four_extincts)-1)/dt, four_extincts[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated Five-agent")
ax2.plot(np.linspace(0, tmax, tmax), four_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments Five-agent')
ax2.set_xlabel('Time')
ax2.set_ylabel('Probability')
ax1.set_ylabel('Probability')

# three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
# five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
ax1.legend(handles=[sim_line, mom_line])
ax1.set_title("Three-agent")
ax2.set_title("Five-agent")
# plt.tight_layout()plt.show()
# plt.show()
plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/prob_of_extinction_plot_both_layers-sims'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
plt.close()


############################## Prob extinct - same plot ##############################################

temp_t = 6000
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
    four_extinction_mom_b_withdead_long = np.load(f)
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
    two_extinction_mom_b_withdead_long = np.load(f)

fig, axs = plt.subplots(2, 2, sharey = True)

axs[0,0].plot(np.arange(0,len(two_extincts)-1)/dt, two_extincts[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated")
axs[0,0].plot(np.linspace(0, tmax, tmax), two_extinction_mom_b_withdead, c = CB_color_cycle[0], linewidth = 2, label = 'Method of Moments')
axs[1,0].plot(np.arange(0,len(four_extincts)-1)/dt, four_extincts[:-1]/num_sims, c = CB_color_cycle[4], ls = 'dotted', linewidth = 3, label= "Simulated")
axs[1,0].plot(np.linspace(0, tmax, tmax), four_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
axs[1,0].set_xlabel('Time')
axs[1,0].set_ylabel('Probability')


axs[0,1].plot(np.linspace(0, temp_t, temp_t), two_extinction_mom_b_withdead_long, c = CB_color_cycle[0], linewidth = 2 )
axs[1,1].plot(np.linspace(0, temp_t, temp_t), four_extinction_mom_b_withdead_long, c = CB_color_cycle[4], linewidth = 2)
axs[1,1].set_xlabel('Time')
axs[0,0].set_ylabel('Probability')

mom_line = Line2D([0], [0], color="black", linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color="black", linewidth=2, linestyle='dotted', label = "Simulations")
three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent') 
fig.legend(handles=[sim_line, mom_line, three_patch, five_patch], ncol = 4)

fig.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/prob_of_extinction_plot_both_layers_2by2-sims'+str(num_sims)+'-mixtimes_'+date+'.pdf', format = 'pdf')


 


############################### AVERAGE BASAL HISTORY - TWO ###########################################
date = datetime.today().strftime('%m-%d-%Y')

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.plot(np.arange(0,len(two_basal_history_other)-1)/dt, two_basal_history_other[1:]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax1.plot(np.linspace(0,tmax, num= tmax), two_basal_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
ax1.set_ylabel("Average Basal")

############################### AVERAGE PARABASAL HISTORY ###########################################


ax2.plot(np.arange(0,len(two_parabasal_history_other))/dt, two_parabasal_history_other/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax2.plot(np.linspace(0,tmax, tmax), two_para_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
ax2.set_ylabel("Average Parabasal")


############################### AVERAGE DEAD HISTORY ###########################################


ax3.plot(np.arange(0,len(two_dead_count))/dt, two_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax3.plot(np.linspace(0,tmax, tmax), two_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
ax3.set_ylabel("Average Dead")

mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
ax1.legend(handles=[sim_line, mom_line])
ax1.set_title("Three-agent Averages")
# plt.show()
fig.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/avg_basal_and_parabasal_dead_plot-sims'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
plt.close()


############################### AVERAGE BASAL HISTORY - FOUR ###########################################
date = datetime.today().strftime('%m-%d-%Y')

fig, axs = plt.subplots(5, sharex=True, figsize=(5, 7.5))

axs[0].plot(np.arange(0,len(four_basal_history_other)-1)/dt, four_basal_history_other[1:]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[0].plot(np.linspace(0,tmax, num= tmax), four_basal_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[0].set_ylabel("Average Basal")

############################### AVERAGE PARABASAL HISTORY - FOUR ###########################################


axs[1].plot(np.arange(0,len(four_parabasal_history_other))/dt, four_parabasal_history_other/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[1].plot(np.linspace(0,tmax, tmax), four_para_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[1].set_ylabel("Average Parabasal")

############################### AVERAGE INTERMEDIATE HISTORY - FOUR ###########################################


axs[2].plot(np.arange(0,len(four_intermed_history))/dt, four_intermed_history/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[2].plot(np.linspace(0,tmax, tmax), four_intermed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[2].set_ylabel("Average Intermediate")

############################### AVERAGE SUPER HISTORY - FOUR ###########################################


axs[3].plot(np.arange(0,len(four_super_history))/dt, four_super_history/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[3].plot(np.linspace(0,tmax, tmax),four_super_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[3].set_ylabel("Average Super")


############################### AVERAGE DEAD HISTORY - FOUR ###########################################


axs[4].plot(np.arange(0,len(four_dead_count))/dt, four_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[4].plot(np.linspace(0,tmax, tmax), four_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
axs[4].set_ylabel("Average Dead")
mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
axs[0].legend(handles=[sim_line, mom_line])
axs[0].set_title("Five-agent Averages")
# plt.show()
#fig.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/figures/mom-sims/avg_5state_plot-sims'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
#plt.close()

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


############################### AVERAGE VIRIONS - Both Layers ################################################


plt.rcParams["font.family"] = "Times New Roman"

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(np.arange(0,len(two_dead_count))/dt, 1000*two_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
ax1.plot(np.linspace(0, tmax, tmax), 1000*two_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
ax2.plot(np.arange(0,len(four_dead_count))/dt, 1000*four_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
ax2.plot(np.linspace(0, tmax, tmax), 1000*four_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')

ax2.set_xlabel('Time')
ax2.set_ylabel('Average Cumulative Virions')
ax1.set_ylabel('Average Cumulative Virions')


# three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
# five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
ax1.legend(handles=[sim_line, mom_line])
ax1.set_title("Three-agent")
ax2.set_title("Five-agent")
# plt.show()
plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/avg_virons_plot-sims_both_layer'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
plt.close()




############################### AVERAGE VIRIONS - both on same plot ################################################

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
fig, ax = plt.subplots()
axins = inset_axes(ax, loc=2, width="40%", height="40%")
axins.tick_params(labelleft=False, labelbottom=False)

three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
mom_line = Line2D([0], [0], color="black", linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color="black", linewidth=2, linestyle='dotted', label = "Simulations") 
ax.legend(handles=[three_patch, five_patch, sim_line, mom_line], loc = 4)

import scipy.stats as stats
diff_dead = 1000*(four_dead_count-two_dead_count)/num_sims
slope, intercept, r, p, std_err = stats.linregress(np.arange(0,len(diff_dead))/dt, diff_dead)
x = np.arange(0,len(diff_dead))/dt
y = diff_dead
line = slope * x + intercept
axins.plot(x, y, alpha = 0.95, c = CB_color_cycle[5])
axins.plot(x, line, c = 'grey')
axins.fill_between(x, line+1.96*std_err, line-1.96*std_err)


ax.plot(np.arange(0,len(two_dead_count))/dt, 1000*two_dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "three Simulated - %d simulations" %(num_sims))
ax.plot(np.linspace(0, tmax, tmax), 1000*two_shed_first_momentwithdead, c = CB_color_cycle[0], linewidth = 2, label = ' three Method of Moments')
ax.plot(np.arange(0,len(four_dead_count))/dt, 1000*four_dead_count/num_sims, c = CB_color_cycle[4], ls = 'dotted', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
ax.plot(np.linspace(0, tmax, tmax), 1000*four_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')

ax.set_xlabel('Time')
ax.set_ylabel('Average Cumulative Virions')
plt.show()




# plt.plot(np.arange(0,len(two_dead_count))/dt, 1000*(four_dead_count - two_dead_count)/num_sims, c = CB_color_cycle[1])
# plt.plot(np.arange(0,len(two_dead_count))/dt,np.arange(0,len(two_dead_count))/dt, c = "black")
# plt.xlim(0, 750)
# plt.xlabel("Time")
# plt.ylabel("Difference in Average Cumulative Virions")





############################### AVERAGE BASAL FOR PERSISTENT INFECTIONS ################################################

# Sims

two_basal_history_nonextinct[two_basal_history_nonextinct == 0] = np.nan
total_non_extincts = np.nanmean(two_basal_history_nonextinct, axis = 0)


# Analytical
mom_nonextincts = np.zeros((tmax))
for ti in range(0, tmax):
    p = ((2*two_basal_first_momentwithdead[ti])/(two_basal_first_momentwithdead[ti] + two_basal_second_momentwithdead[ti]))
    mom_nonextincts[ti] = 1/p


plt.plot(np.arange(0,len(total_non_extincts)-1)/dt, total_non_extincts[1:], c = CB_color_cycle[0], ls = 'dashed', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0,tmax-1, num = tmax-1), mom_nonextincts[1:], c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
#plt.fill_between(np.arange(0,len(total_non_extincts_100)-1)/10, total_non_extincts_100[1:]+ var_total_non_extincts_100[1:], total_non_extincts_100[1:], alpha = 0.5)
plt.legend(loc= 'best')
plt.xlabel('Time')
plt.ylabel('Basals excluding extinction')
plt.tight_layout()
plt.show()
#plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/figures/mom-sims/avg_nonextinct_basals_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()



############################### Double checks on figures ################################################
temp_t = 6000
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
    four_extinction_mom_b_withdead_long = np.load(f)
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
    two_extinction_mom_b_withdead_long = np.load(f)

# Prob of Extinction

plt.plot(np.linspace(0, temp_t, temp_t), two_extinction_mom_b_withdead_long, ls = 'dotted', c = CB_color_cycle[3], linewidth = 3, label = 'Method of Moments Three-agent')
plt.plot(np.linspace(0, temp_t, temp_t), four_extinction_mom_b_withdead_long, c = CB_color_cycle[0], linewidth = 2, label = 'Method of Moments Five-agent')
plt.xlabel('Time')
plt.ylabel('Probability')

plt.legend()
plt.tight_layout()
plt.show()