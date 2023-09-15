#Plotting for Two state 
import numpy as np
import matplotlib.pyplot as plt

with open('cumu_extinct_2layer_main_delta.npy', 'rb') as f:
    cumu_extinct_delta = np.load(f)
with open('cumu_extinct_2layer_main_poisson.npy', 'rb') as f:
    cumu_extinct_poisson= np.load(f)
    
with open('extinct_mom_b_2layer_main_delta', 'rb') as f:
    extinct_mom_b_delta = np.load(f)

with open('extinct_mom_b_2layer_main_poisson', 'rb') as f:
    extinct_mom_b_poisson = np.load(f)


t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)













plt.plot(t_vec, cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta Approx)')
plt.plot(t_vec, extinct_mom_b_delta, label = 'Cumulative probability of extinction - MOM basals (Delta Approx)')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('2 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()


plt.plot(t_vec, cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta Approx)')
plt.plot(t_vec, extinct_mom_b_poisson, label = 'Cumulative probability of extinction - MOM basals (Poisson Approx)')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('2 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()


# import matplotlib.ticker as ticker
# # # # # # Plot
# majors = []
# m_path_vals = []
# y_vals = []
# z_vals = []
# ax = plt.figure().add_subplot(projection='3d')
# stddev = np.sqrt(m_path[-1][2]-m_path[-1][0]**2)
# for t in range(0, t_length-1, 100):
#     ax.plot(range(nb_of_states_b), np.sum(x_path[t], axis=1), t, marker="o", ls='--')
#     majors.append(t)
#     m_path_vals.append(m_path[t][0])
#     y_vals.append(0)
#     z_vals.append(t)
# ax.stem(m_path_vals, y_vals, z_vals, linefmt='black', markerfmt = 'black', label='First moment')
# #ax.stem(m_path_vals+stddev, y_vals, z_vals, linefmt='', label='Standard deviation')
# #ax.legend()
# ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
# ax.set_ylabel('Probability')
# ax.set_xlabel('Number of infected basal cells')
# ax.set_zlabel('Time')
# ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
# plt.show()
# #plt.savefig("explicit_basal_2layer_3D.pdf", format = "pdf")
# plt.close()

# # # #
# # # # # Plot
# ax = plt.figure().add_subplot(projection='3d')
# majors = []
# m_path_vals = []
# y_vals = []
# z_vals = []
# stddev = np.sqrt(m_path[-1][3]-m_path[-1][1]**2)
# for t in range(0,t_length-1,100):
#     ax.plot(range(nb_of_states_p), np.sum(x_path[t], axis=0), t, marker="o", ls='--')
#     majors.append(t)
#     m_path_vals.append(m_path[t][0])
#     y_vals.append(0)
#     z_vals.append(t)
# ax.stem(m_path_vals, y_vals, z_vals, linefmt='black', markerfmt = 'black', label='First moment')
# #ax.stem(m_path_vals+stddev, y_vals, z_vals, linefmt='--', label='Standard deviation')
# ax.view_init(elev=51, azim=39, roll=120) # Elevation is bottom spin, Azimth is twist R-L 
# ax.set_ylabel('Probability')
# ax.set_xlabel('Number of infected parabasal cells')
# ax.set_zlabel('Time')
# ax.zaxis.set_major_locator(ticker.FixedLocator(majors))
# plt.show()
# #plt.savefig("explicit_parabasal_2layer_3D.pdf", format = "pdf")
# plt.close()





# plt.vlines(m_path[1][0], 0, 0.9, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[1][0]+stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# #plt.vlines(m_path[1][0]-stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# plt.plot(range(nb_of_states_b), np.sum(x_path[1], axis = 1), marker="o", lw=0, ms=15, label=fr"$t = {t_vec[1]}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected basal cells')



# plt.savefig("MoM_basal_compare_time1.pdf", format = "pdf")
# plt.close()

# plt.vlines(m_path[1][1], 0, 0.3, colors='black', linestyles='-', label='First moment')
# #plt.vlines(m_path[1][1]+stddev_parabasal, 0, 0.175, colors='gray', linestyles='--', label='Standard deviation')
# #plt.vlines(m_path[1][1]-stddev_parabasal, 0, 0.325, colors='gray', linestyles='--', label='Standard deviation')
# plt.plot(range(nb_of_states_p), np.sum(x_path[1], axis=0), marker="o", lw=0, ms=15, label=fr"$t = {t_vec[1]}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected parabasal cells')

# plt.savefig("MoM_parabasal_compare_time1.pdf", format = "pdf")
# plt.close()
