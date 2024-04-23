# import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math
import csv



###################### PLOTTING ##############################

# Input values
num_sims = 50000
tmax = 500
dt = 100


###################### OPENING SIMULATION FILES ##############################


from datetime import datetime
date = datetime.today().strftime('%m-%d-%Y')
date = "03-15-2024"


# with open('extinction_times_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     extinction_times = np.load(f)

with open('four_layer/data/otherbasal_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    basal_history_other = np.load(f) 

with open('four_layer/data/otherparabasal_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    parabasal_history_other = np.load(f)

with open('four_layer/data/intermed_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    intermed_history = np.load(f)

with open('four_layer/data/super_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    super_history = np.load(f)

with open('four_layer/data/otherdead_history_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    dead_count = np.load(f)

# with open('otherparabasal_variable_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
#     parabasal_variable = np.load(f)


with open('four_layer/data/otherextinct_'+date+'_time500_sims%d.npy'%(num_sims), 'rb') as f:
    extincts = np.load(f)



###################### OPENING MOM FILES ##############################



with open('four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'rb') as handle:    
     extinct_mom_b_geometric = np.load(handle)





# with open('cumu_extinct_2layer_main_delta_1_8.npy', 'rb') as f:
#     cumu_extinct_delta = np.load(f)

# with open('shed_first_moments_delta_02-23_time500.npy', 'rb') as f:
#     shed_first_moment_delta = np.load(f)

# with open('basal_first_moment_geom_1-29_500.npy', 'rb') as f:
#     basal_first_moment= np.load(f)
# with open('para_first_moment_geom_1-29_500.npy', 'rb') as f:
#     para_first_moment = np.load(f)
# with open('para_second_moment_geom_02-29_time500.npy', 'rb') as f:
#     para_second_moment = np.load(f)


with open('four_layer/data/extinction_mom_b_4layerdead_geometric_'+date+'_time500.npy', 'rb') as f:
    extinction_mom_b_withdead = np.load(f)

with open('four_layer/data/shed_first_moments_delta_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    shed_first_momentwithdead = np.load(f)

with open('four_layer/data/basal_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    basal_first_momentwithdead = np.load(f)

with open('four_layer/data/para_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    para_first_momentwithdead = np.load(f)

with open('four_layer/data/intermed_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    intermed_first_momentwithdead = np.load(f)

with open('four_layer/data/super_first_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    super_first_momentwithdead = np.load(f)

with open('four_layer/data/para_second_moment_geom_4layerwithdead_'+date+'_time500.npy', 'rb') as f:
    para_second_momentwithdead = np.load(f)


############################### COLORS ##################################
    

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


############################### PROBABILITY OF EXTINCTION ###########################################



plt.plot(np.arange(0,len(extincts)-1)/dt, extincts[:-1]/num_sims, c = CB_color_cycle[0], ls = 'dashdot', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0, tmax, tmax), extinction_mom_b_withdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
plt.title("Probability of Extinction of Infected Basal Cells")
plt.xlabel('Time')
plt.ylabel('Probability')
plt.legend()
#plt.show()
plt.savefig('four_layer/figures/mom-sims/prob_of_extinction_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()


############################### AVERAGE BASAL HISTORY ###########################################
date = datetime.today().strftime('%m-%d-%Y')

fig, axs = plt.subplots(5, sharex=True)

axs[0].plot(np.arange(0,len(basal_history_other)-1)/dt, basal_history_other[1:]/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[0].plot(np.linspace(0,tmax, num= tmax), basal_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[0].legend(loc= 'best')
axs[0].set_ylabel("Basal")

############################### AVERAGE PARABASAL HISTORY ###########################################


axs[1].plot(np.arange(0,len(parabasal_history_other))/dt, parabasal_history_other/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[1].plot(np.linspace(0,tmax, tmax), para_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[1].set_ylabel("Parabasal")

############################### AVERAGE INTERMEDIATE HISTORY ###########################################


axs[2].plot(np.arange(0,len(intermed_history))/dt, intermed_history/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[2].plot(np.linspace(0,tmax, tmax), intermed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[2].set_ylabel("Intermediate")

############################### AVERAGE SUPER HISTORY ###########################################


axs[3].plot(np.arange(0,len(super_history))/dt, super_history/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[3].plot(np.linspace(0,tmax, tmax), super_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
axs[3].set_ylabel("Super")


############################### AVERAGE DEAD HISTORY ###########################################


axs[4].plot(np.arange(0,len(dead_count))/dt, dead_count/num_sims, c = CB_color_cycle[0], ls = 'dotted', linewidth = 3, label = "Simulated - %d simulations" %(num_sims))
axs[4].plot(np.linspace(0,tmax, tmax), shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = "Method of Moments")
plt.xlabel("Time")
axs[4].set_ylabel("Dead")

axs[0].set_title("Averages of States")
#plt.show()
#fig.savefig('four_layer/figures/mom-sims/avg_5state_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
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


############################### AVERAGE VIRIONS ################################################

plt.plot(np.arange(0,len(dead_count))/dt, 1000*dead_count/num_sims, c = CB_color_cycle[0], ls = 'dashdot', linewidth = 3, label= "Simulated - %d simulations" %(num_sims))
plt.plot(np.linspace(0, tmax, tmax), 1000*shed_first_momentwithdead, c = CB_color_cycle[4], linewidth = 2, label = 'Method of Moments')
plt.title('Average Virions over Time')
plt.xlabel('Time')
plt.ylabel('Average Virons')
plt.legend()
#plt.show()
plt.savefig('four_layer/figures/mom-sims/avg_virons_plot-sims'+str(num_sims)+'-time500_'+date+'.pdf', format = 'pdf')
plt.close()

