#Plotting for Four state 
import numpy as np
import matplotlib.pyplot as plt

with open('x_path_file_4layer_main.npy', 'rb') as f:
    x_path = np.load(f)
    
with open('m_path_4layer_main_1.npy', 'rb') as f:
    m_path = np.load(f)

t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)

prop = 0.3

nb_of_states_b = int(100*prop)
nb_of_states_p = int(500*prop)
nb_of_states_i = int(50*prop)
nb_of_states_s = int(50*prop)


## Solving for extinction probabilities 
# ax.plot(range(nb_of_states_b), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 1), axis = 1), t, marker="o", ls='--')

## Explicit ME

extinct_b = np.zeros((t_length))
cumu_extinct = np.zeros((t_length))
for t in range(t_length - 1):
    t_plus = np.sum(np.sum(np.sum(x_path[t+1], axis = 0), axis = 1), axis = 1)
    t_current = np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 1), axis = 1)
    extinct_b[t] = t_plus - t_current ## TODO put in correct summations 
    cumu_extinct[t] = np.sum(extinct_b) #x_path[t][0][0]

# MOM


extinct_mom_b = np.zeros((t_length))
for t in range(t_length - 1):
    extinct_mom_b[t] = (1-((m_path[t][0]**2)/(m_path[t][4])))


plt.plot(t_vec, cumu_extinct, label = 'Cumulative probability of extinction - Explicit')
plt.plot(t_vec, extinct_mom_b, label = 'Cumulative probabilty of extinction - MOM basals')
# plt.plot(t_vec, extinct_mom_b, label = 'Cumulative probability of extinction - MOM both ')
plt.legend()
plt.ylabel('Probability')
plt.xlabel('Time')
plt.title('4 Layer system')
# plt.savefig("extinction_prob.pdf", format = 'pdf')
plt.show()
plt.close()




# # # # # Plot
# for t in np.arange(1, t_length-1, 100):
#     plt.plot(range(nb_of_states_b), np.sum(np.sum(np.sum(x_path[t], axis = 1), axis = 1), axis = 1), marker="o", ls='--',
#              label=fr"$t = {t_length * t / t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected basal cells')
# plt.show()
# #plt.savefig("flat_explicit_basal_4Layer.pdf", format = "pdf")
# plt.close()
# # # #
# # # # # Plot
# for t in np.arange(1,t_length-1,100):
#     plt.plot(range(nb_of_states_p), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 1), axis = 1), marker="o", ls='--', label=fr"$t = {t_length*t/t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected parabasal cells')
# plt.show()
# #plt.savefig("explicit_parabasal.pdf", format = "pdf")
# plt.close()

# # # # # Plot
# for t in np.arange(1,t_length-1,100):
#     plt.plot(range(nb_of_states_i), np.sum(np.sum(np.sum(x_path[t], axis = 0), axis = 0), axis = 1), marker="o", ls='--', label=fr"$t = {t_length*t/t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected intermediate cells')
# plt.show()
# #plt.savefig("explicit_parabasal.pdf", format = "pdf")
# plt.close()

# # # # # Plot
# for t in np.arange(1,t_length-1,100):
#     plt.plot(range(nb_of_states_s), np.sum(np.sum(np.sum(x_path[t], axis = 1), axis = 1), axis = 1), marker="o", ls='--', label=fr"$t = {t_length*t/t_steps:.2f}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected superficial cells')
# plt.show()
# #plt.savefig("explicit_parabasal.pdf", format = "pdf")
# plt.close()

# stddev = np.sqrt(m_path[-1][4]-m_path[-1][0]**2)
# stddev_parabasal = np.sqrt(m_path[-1][5]-m_path[-1][1]**2)


# plt.vlines(m_path[10][0], 0, 0.9, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[10][0]+stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# #plt.vlines(m_path[1][0]-stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# plt.plot(range(nb_of_states_b), np.sum(np.sum(np.sum(x_path[10], axis = 1), axis = 1), axis = 1), marker="o", lw=0, ms=15, label=fr"$t = {t_vec[10]}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected basal cells')
# plt.show()
# plt.close()

# plt.vlines(m_path[10][1], 0, 0.9, colors='black', linestyles='-', label='First moment')
# plt.vlines(m_path[10][1]+stddev_parabasal, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# #plt.vlines(m_path[1][0]-stddev, 0, 0.9, colors='gray', linestyles='--', label='Standard deviation')
# plt.plot(range(nb_of_states_p), np.sum(np.sum(np.sum(x_path[10], axis = 0), axis = 1), axis = 1), marker="o", lw=0, ms=15, label=fr"$t = {t_vec[10]}$")
# plt.legend()
# plt.ylabel('Occupation number')
# plt.xlabel('Number of infected parabasal cells')
# plt.show()
# plt.close()

# plt.plot(t_vec, extinct, label = 'Probability of extinction')
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()

# plt.plot(t_vec, cumu_extinct, label = 'Cumulative probability of extinction')
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()