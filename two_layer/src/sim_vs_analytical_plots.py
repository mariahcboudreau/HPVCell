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

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/otherbasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_basal_history_other_100 = np.load(f) 

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/total_nonextincts_05-04-2024_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_total_nonextincts = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/two_dead_history_mean_05-06-2024_time%d_sims%d.npy'%(tmax, num_sims), 'rb') as f:
    two_dead_history_mean = np.load(f)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/otherparabasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_parabasal_history_other_100 = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/otherdead_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_dead_count_100 = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/otherextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    two_extincts_100 = np.load(f)


##### Four Layer data 

# with open('extinction_times_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     extinction_times = np.load(f)

date = "04-23-2024"

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/otherbasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_basal_history_other_100 = np.load(f) 

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/otherparabasal_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_parabasal_history_other_100  = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/intermed_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_intermed_history_100  = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/super_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_super_history_100  = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/otherdead_history_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_dead_count_100  = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/otherextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_extincts_100  = np.load(f)

date = "05-06-2024"


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/total_nonextincts_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'rb') as f:
    four_total_non_extincts_10 = np.load(f)

date = "05-13-2024"
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/diff_dead_history_mean_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'rb') as f:
    diff_dead_mean = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/diff_dead_history_variance_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'rb') as f:
    diff_dead_var = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/four_dead_history_nonextinct_mean_'+date+'_time%d_sims%d.npy'%(tmax, num_sims), 'rb') as f:
    four_dead_history_nonextinct_mean_10 = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/two_dead_history_nonextinct_mean_05-13-2024_time%d_sims%d.npy'%(tmax, num_sims), 'rb') as f:
    two_dead_history_nonextinct_mean = np.load(f)


###################### OPENING MOM FILES ##############################

date = "05-13-2024"

# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'rb') as handle:    
#      extinct_mom_b_geometric = np.load(handle)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_extinction_mom_b_withdead = np.load(f)


with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/shed_first_moments_delta_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_shed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/shed_second_moments_delta_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_shed_second_momentwithdead = np.load(f)  

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/basal_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_basal_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/para_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_para_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/intermed_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_intermed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/super_first_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_super_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/basal_second_moment_geom_4layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    four_basal_second_momentwithdead = np.load(f)




# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCellSim/two_layer/data/extinction_mom_b_2layer_geometric_'+date+'_time500.npy', 'rb') as handle:    
#      extinct_mom_b_geometric = np.load(handle)



with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_extinction_mom_b_withdead = np.load(f)



with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/shed_first_moments_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_shed_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/shed_second_moments_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_shed_second_momentwithdead = np.load(f)  

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/basal_first_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_basal_first_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/basal_second_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_basal_second_momentwithdead = np.load(f)

with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/para_first_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
    two_para_first_momentwithdead = np.load(f)

# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/para_second_moment_geom_2layerwithdead_'+date+'_time%d.npy'%(tmax), 'rb') as f:
#     two_para_second_momentwithdead = np.load(f)


temp_t = 6000
date = '04-25-2024'
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
    four_extinction_mom_b_withdead_long = np.load(f)
with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
    two_extinction_mom_b_withdead_long = np.load(f)

############################### COLORS/IMPORTS ##################################
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
date = datetime.today().strftime('%m-%d-%Y')
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["figure.dpi"] = 200

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


############################### PROBABILITY OF EXTINCTION - Both Layers ###########################################


fig, (ax1, ax2) = plt.subplots(2, sharex=True)


ax1.plot(np.arange(0,len(two_extincts_100)-1)/dt, two_extincts_100[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated Three-agent")
ax1.plot(np.linspace(0, tmax, tmax), two_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments Three-agent')
ax2.plot(np.arange(0,len(four_extincts_100)-1)/dt, four_extincts_100[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated Five-agent")
ax2.plot(np.linspace(0, tmax, tmax), four_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments Five-agent')
ax2.set_xlabel('Time')
ax2.set_ylabel('Probability')
ax1.set_ylabel('Probability')

# three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
# five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
ax1.legend(handles=[sim_line, mom_line], frameon = False)
ax1.set_title("Three-agent")
ax2.set_title("Five-agent")
# plt.tight_layout()plt.show()
plt.show()
#plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/prob_of_extinction_plot_both_layers-sims'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
plt.close()


############################## Prob extinct - same plot ##############################################



fig, (ax1, ax2) = plt.subplots(1,2, sharey = True, figsize = (10,4))
axins = inset_axes(ax2, loc=4, width="40%", height="40%")
axins.tick_params(labelleft=False, labelbottom=False)

x = np.linspace(0, temp_t, temp_t)
y1 = two_extinction_mom_b_withdead_long
y2 = four_extinction_mom_b_withdead_long
slope_three = []
for i in range(2, len(x),2):
    slope_three.append((y1[i]-y1[i-1])/(x[i]-x[i-1]))


ax1.plot(np.arange(0,len(two_extincts_100)-1)/dt, two_extincts_100[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 6, alpha = 0.7, label= "Simulated")
ax1.plot(np.linspace(0, tmax, tmax), two_extinction_mom_b_withdead, c = CB_color_cycle[0], linewidth = 6, alpha = 0.5, label = 'Method of Moments')
ax1.plot(np.arange(0,len(four_extincts_100)-1)/dt, four_extincts_100[:-1]/num_sims, c = CB_color_cycle[4], ls = 'dotted', linewidth = 3, label= "Simulated")
ax1.plot(np.linspace(0, tmax, tmax), four_extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
ax1.set_xlabel('Time')
ax1.set_ylabel('Probability')


ax2.plot(np.linspace(0, temp_t, temp_t), two_extinction_mom_b_withdead_long, c = CB_color_cycle[0], linewidth = 6, alpha = 0.5)
ax2.plot(np.linspace(0, temp_t, temp_t), four_extinction_mom_b_withdead_long, c = CB_color_cycle[4], linewidth = 2)
ax2.set_xlabel('Time')

mom_line = Line2D([0], [0], color="black", linewidth=3, label = "Method of moments") 
sim_line = Line2D([0], [0], color="black", linewidth=3, linestyle='dotted', label = "Simulations")
three_patch = mpatches.Patch(color=CB_color_cycle[0], alpha = 0.5, label='Three-agent')
five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent') 
fig.legend(handles=[sim_line, mom_line, three_patch, five_patch],loc = 9, ncol = 4, frameon = False)

fig.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/prob_of_extinction_plot_both_layers_2by2-sims'+str(num_sims)+'-mixtimes_'+date+'.pdf', format = 'pdf')


 


############################### AVERAGE BASAL HISTORY - TWO ###########################################


fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.plot(np.arange(0,len(two_basal_history_other_100)-1)/dt, two_basal_history_other_100[1:]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax1.plot(np.linspace(0,tmax, num= tmax), two_basal_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
ax1.set_ylabel("Avg. Basal")

############################### AVERAGE PARABASAL HISTORY ###########################################


ax2.plot(np.arange(0,len(two_parabasal_history_other_100))/dt, two_parabasal_history_other_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax2.plot(np.linspace(0,tmax, tmax), two_para_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
ax2.set_ylabel("Avg. Parabasal")


############################### AVERAGE DEAD HISTORY ###########################################


ax3.plot(np.arange(0,len(two_dead_count_100))/dt, two_dead_count_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
ax3.plot(np.linspace(0,tmax, tmax), two_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
ax3.set_ylabel("Avg. Dead")

mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
ax1.legend(handles=[sim_line, mom_line], frameon = False)
# plt.show()
fig.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/avg_basal_and_parabasal_dead_plot-sims'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
plt.close()


############################### AVERAGE BASAL HISTORY - FOUR ###########################################

fig, axs = plt.subplots(5, sharex=True, figsize = (6.4, 6.4))

axs[0].plot(np.arange(0,len(four_basal_history_other_100)-1)/dt, four_basal_history_other_100[1:]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[0].plot(np.linspace(0,tmax, num= tmax), four_basal_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[0].set_ylabel("Avg. Basal")

############################### AVERAGE PARABASAL HISTORY - FOUR ###########################################


axs[1].plot(np.arange(0,len(four_parabasal_history_other_100))/dt, four_parabasal_history_other_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[1].plot(np.linspace(0,tmax, tmax), four_para_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[1].set_ylabel("Avg. Parabasal")

############################### AVERAGE INTERMEDIATE HISTORY - FOUR ###########################################


axs[2].plot(np.arange(0,len(four_intermed_history_100))/dt, four_intermed_history_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[2].plot(np.linspace(0,tmax, tmax), four_intermed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[2].set_ylabel("Avg. Intermediate")

############################### AVERAGE SUPER HISTORY - FOUR ###########################################


axs[3].plot(np.arange(0,len(four_super_history_100))/dt, four_super_history_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[3].plot(np.linspace(0,tmax, tmax),four_super_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[3].set_ylabel("Avg. Super")


############################### AVERAGE DEAD HISTORY - FOUR ###########################################


axs[4].plot(np.arange(0,len(four_dead_count_100))/dt, four_dead_count_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[4].plot(np.linspace(0,tmax, tmax), four_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
axs[4].set_ylabel("Avg. Dead")
mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
fig.legend(handles=[sim_line, mom_line], loc= 9,ncol = 2, frameon = False)
# axs[0].set_title("Five-agent Averages")
# plt.show()
fig.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/figures/mom-sims/avg_5state_plot-sims'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
plt.close()

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


# plt.rcParams["font.family"] = "Times New Roman"

# fig, (ax1, ax2) = plt.subplots(2, sharex=True)
# ax1.plot(np.arange(0,len(two_dead_count_100))/dt, 1000*two_dead_count_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
# ax1.plot(np.linspace(0, tmax, tmax), 1000*two_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
# ax2.plot(np.arange(0,len(four_dead_count_100))/dt, 1000*four_dead_count_100/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
# ax2.plot(np.linspace(0, tmax, tmax), 1000*four_shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')

# ax2.set_xlabel('Time')
# ax2.set_ylabel('Average Cumulative Virions')
# ax1.set_ylabel('Average Cumulative Virions')


# # three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
# # five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
# mom_line = Line2D([0], [0], color=CB_color_cycle[4], linewidth=2, label = "Method of moments") 
# sim_line = Line2D([0], [0], color=CB_color_cycle[0], linewidth=2, linestyle='dotted', label = "Simulations") 
# ax1.legend(handles=[sim_line, mom_line], frameon = False)
# ax1.set_title("Three-agent")
# ax2.set_title("Five-agent")
# # plt.show()
# plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/avg_virons_plot-sims_both_layer'+str(num_sims)+'-time'+str(tmax)+'_'+date+'.pdf', format = 'pdf')
# plt.close()




############################### AVERAGE VIRIONS - inset plot ################################################


fig, ax = plt.subplots()
axins = inset_axes(ax, loc=2, width="40%", height="40%")
axins.tick_params(labelleft=False, labelbottom=False)


three_patch = mpatches.Patch(color=CB_color_cycle[0], label='Three-agent')
five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
mom_line = Line2D([0], [0], color="black", linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color="black", linewidth=2, linestyle='dotted', label = "Simulations") 
diff_line = Line2D([0], [0], color=CB_color_cycle[5], linewidth=2, label = "Difference between systems") 
slope_line = Line2D([0], [0], color='grey', linewidth=2, label = "Linear Regression") 
fig.legend(handles=[diff_line, slope_line, three_patch, five_patch, sim_line, mom_line], ncol = 3, loc = 9, frameon = False)

import scipy.stats as stats
slope, intercept, r, p, std_err = stats.linregress(np.arange(0,len(diff_dead_mean))/10, diff_dead_mean)
x = np.arange(0,len(diff_dead_mean))/10
y = diff_dead_mean
line = slope * x + intercept
temp = y+ 0.5*diff_dead_var
# axins.plot(x,temp, alpha = 0.5, c = CB_color_cycle[5])
axins.plot(x, y, alpha = 0.95, c = CB_color_cycle[5])
axins.plot(x, line, c = 'grey')
# axins.fill_between(x, y, temp, alpha = 0.5, color = CB_color_cycle[5])

# axins.fill_between(x, line+1.96*std_err, line-1.96*std_err)


## MoM estimated non-extinct values 
# Analytical
two_mom_shed_nonextincts = np.zeros((tmax))
for ti in range(0, tmax):
    p = ((2*two_shed_first_momentwithdead[ti])/(two_shed_first_momentwithdead[ti] + two_shed_second_momentwithdead[ti]))
    two_mom_shed_nonextincts[ti] = 1/p

# Analytical
four_mom_shed_nonextincts = np.zeros((tmax))
for ti in range(0, tmax):
    p = ((2*four_shed_first_momentwithdead[ti])/(four_shed_first_momentwithdead[ti] + four_shed_second_momentwithdead[ti]))
    four_mom_shed_nonextincts[ti] = 1/p


ax.plot(np.arange(0,len(two_dead_history_nonextinct_mean))/10, 1000*two_dead_history_nonextinct_mean, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label= "three Simulated - %d simulations" %(num_sims))
ax.plot(np.linspace(0, tmax, tmax), 1000*two_shed_first_momentwithdead, c = CB_color_cycle[0], linewidth = 2, label = ' three Method of Moments')
ax.plot(np.arange(0,len(four_dead_history_nonextinct_mean_10))/10, 1000*four_dead_history_nonextinct_mean_10, c = CB_color_cycle[4], ls = 'dotted', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
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

# two_basal_history_nonextinct_10[two_basal_history_nonextinct_10 == 0] = np.nan
# total_non_extincts = np.nanmean(two_basal_history_nonextinct_10, axis = 0)


date = datetime.today().strftime('%m-%d-%Y')
# Analytical
two_mom_nonextincts = np.zeros((tmax))
for ti in range(0, tmax):
    p = ((2*two_basal_first_momentwithdead[ti])/(two_basal_first_momentwithdead[ti] + two_basal_second_momentwithdead[ti]))
    two_mom_nonextincts[ti] = 1/p

# Analytical
four_mom_nonextincts = np.zeros((tmax))
for ti in range(0, tmax):
    p = ((2*four_basal_first_momentwithdead[ti])/(four_basal_first_momentwithdead[ti] + four_basal_second_momentwithdead[ti]))
    four_mom_nonextincts[ti] = 1/p


plt.plot(np.arange(0,len(two_total_nonextincts)-1)/10, two_total_nonextincts[1:], c = CB_color_cycle[0], alpha = 0.4, ls = 'dotted', linewidth = 6, label = "Simulations")
plt.plot(np.linspace(0,tmax-1, num = tmax-1), two_mom_nonextincts[1:], c = CB_color_cycle[0], linewidth = 5, alpha = 0.5, label = "Method of Moments")
plt.plot(np.arange(0,len(four_total_non_extincts_10)-1)/10, four_total_non_extincts_10[1:], c = CB_color_cycle[4], ls = 'dotted', linewidth = 3, label = "Simulations")
plt.plot(np.linspace(0,tmax-1, num = tmax-1), four_mom_nonextincts[1:], c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")

three_patch = mpatches.Patch(color=CB_color_cycle[0], alpha = 0.6, label='Three-agent')
five_patch = mpatches.Patch(color=CB_color_cycle[4], label='Five-agent')
mom_line = Line2D([0], [0], color="black", linewidth=2, label = "Method of moments") 
sim_line = Line2D([0], [0], color="black", linewidth=2, linestyle='dotted', label = "Simulations") 
plt.legend(handles=[three_patch, five_patch, sim_line, mom_line], ncol = 2, loc = 2, frameon = False)

plt.xlabel('Time')
plt.ylabel('Average Non-extinct Basals')
plt.tight_layout()
plt.show()
plt.savefig('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/figures/mom-sims/avg_nonextinct_basals_plot-sims-bothagents'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()



############################### Double checks on figures ################################################
# temp_t = 6000
# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
#     four_extinction_mom_b_withdead_long = np.load(f)
# with open('/Users/mcboudre/Documents/MOCS2/testCode/HPVCell/two_layer/data/extinction_mom_b_2layerdead_geometric_'+date+'_time%d.npy'%(temp_t), 'rb') as f:
#     two_extinction_mom_b_withdead_long = np.load(f)

# # Prob of Extinction - long

# plt.plot(np.linspace(0, temp_t, temp_t), two_extinction_mom_b_withdead_long, c = CB_color_cycle[3], linewidth = 3, label = 'Method of Moments Three-agent')
# plt.plot(np.linspace(0, temp_t, temp_t), four_extinction_mom_b_withdead_long, c = CB_color_cycle[3], linewidth = 2, label = 'Method of Moments Five-agent')
# plt.xlabel('Time')
# plt.ylabel('Probability')

# plt.legend(frameon = False)
# plt.tight_layout()
# plt.show()