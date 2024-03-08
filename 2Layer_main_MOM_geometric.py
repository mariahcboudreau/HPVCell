# HPVCellSim Master Equations MOM 4 layer

import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import sys
import scipy



# WE will use the odeint routine
from scipy.integrate import odeint
# With a wrapper to facilitate 2d arrays
from odeintw import odeintw

# MOM

# m[0] first moment for b
# m[1] first moment for p
# m[2] second moment for b
# m[3] second moment for p
# m[4] interaction term bp


#MoM
def MOM(m, t, beta, gamma, delta, rho, theta): #figure out the right way to do the average of separate variables
    dx = 0*m
    # First Moment for b
    dx[0] = (beta - delta)* m[0]
    # First Moment for p
    dx[1] = (rho - theta) * m[1] + (2 * delta + gamma) * m[0]
    # Second Moment for b
    dx[2] = 2*(beta - delta) * m[2] + (beta + delta) * m[0]
    # Second Moment for p
    dx[3] = (theta + rho) * m[1] + (2 * rho - 2 * theta) * m[3] + (2 * gamma + 4 * delta) * m[4] + (gamma + 4 * delta) * m[0]
    # Covariance for bp
    dx[4] = beta * m[4] - theta * m[4] + rho * m[4] + gamma * m[2] - delta * m[4] + 2 * delta * m[2] - 2 * delta * m[0]
    # Shed first moment
    dx[5] = 1000*theta*m[1]
    # Shed second moment
    dx[6] = theta*1000**2*m[3]
    return dx


# Time of observations
t_length = 500
t_steps = 500
t_vec = np.linspace(0, t_length, t_steps)

# Initial conditions

m_0_delta_geo = np.zeros(7)
m_0_delta_geo[0] = 1
m_0_delta_geo[2] = 1

m_0_poisson = np.zeros(7)
m_0_poisson[0] = 1
m_0_poisson[2] = 2




# # Parameters of the model OLD
# symm_div = 0.04
# asymm_div = 1 - 2*symm_div # Asymmetric divisions are more probable, 84% chance
# shed = 2.0*(0.0082) # similar to division but just a little bit different according to the plos epithelial strat paper

# rho = symm_div
# gamma = asymm_div
# delta = symm_div
# beta = symm_div + 0.01
# theta = shed

# Updated parameters
R_b = 0.03          # Division rate of basal - Murall citation
R_p = 0.39          # Division rate of parabasal - Murall citation
symm_div = 0.08     # Symmetric division rate - Clayton
asymm_div = 0.84    # Asymmetric division rate - Clayton

beta = R_b * symm_div  + 0.001      # bbb 
gamma = R_b * asymm_div             # bbp
delta = R_b * symm_div              # bpp
rho = R_p * symm_div                # ppp
theta = 0.67                        # Shed - Murall

# Integration
M = lambda m, t: MOM(m, t, beta, gamma, delta, rho, theta)
m_path_poisson = odeintw(M, m_0_poisson, t_vec)
m_path_delta_geo = odeintw(M, m_0_delta_geo, t_vec)


print("done integrating")


# with open('cumu_extinct_2layer_main_delta_last_12-11.npy', 'rb') as f:
#     cumu_extinct_delta_last = np.load(f)
# with open('cumu_extinct_2layer_main_poisson_last_12-11.npy', 'rb') as f:
#     cumu_extinct_poisson_last = np.load(f)
# with open('cumu_extinct_2layer_main_delta_11_27.npy', 'rb') as f:
#     cumu_extinct_delta = np.load(f)
# with open('cumu_extinct_2layer_main_poisson_11_13.npy', 'rb') as f:
#     cumu_extinct_poisson= np.load(f)



# extinct_mom_b_delta = np.zeros((len(t_vec)))
extinct_mom_b_poisson = np.zeros((len(t_vec)))
extinct_mom_p_poisson = np.zeros((len(t_vec)))
extinct_mom_b_geometric = np.zeros((len(t_vec)))
extinct_mom_p_geometric = np.zeros((len(t_vec)))
# extinct_mom_b_expon = np.zeros((len(t_vec)))
# extinct_mom_b_power = np.zeros((len(t_vec)))
shed_moments = np.zeros((2,len(t_vec))) # First moment is index 0 and second moment is index 1
   # # MOM
for ti in range(0, len(t_vec)):

    # extinct_mom_b_delta[ti] = (1-((m_path_delta_geo[ti][0]**2)/(m_path_delta_geo[ti][2])))
    
    extinct_mom_b_poisson[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][2] - m_path_poisson[ti][0])  + scipy.stats.poisson.pmf(k = 0, mu = ((m_path_poisson[ti][2])/(m_path_poisson[ti][0])-1))
    extinct_mom_p_poisson[ti] = 1 - (m_path_poisson[ti][1]**2)/(m_path_poisson[ti][3] - m_path_poisson[ti][1])  + scipy.stats.poisson.pmf(k = 0, mu = ((m_path_poisson[ti][3])/(m_path_poisson[ti][1])-1))
    extinct_mom_b_geometric[ti] = 1 - (2*m_path_delta_geo[ti][0]**2) / (m_path_delta_geo[ti][2] + m_path_delta_geo[ti][0] )  + scipy.stats.geom.pmf(k = 0, p = (2/(m_path_delta_geo[ti][0]*m_path_delta_geo[ti][2] + 1)))
    extinct_mom_p_geometric[ti] = 1 - (2*m_path_delta_geo[ti][1]**2) / (m_path_delta_geo[ti][3] + m_path_delta_geo[ti][1] )  + scipy.stats.geom.pmf(k = 0, p = (2/(m_path_delta_geo[ti][1]*m_path_delta_geo[ti][3] + 1)))
    # shed_moments[0][ti] = 1000*theta*m_path_delta_geo[ti][1]
    # shed_moments[1][ti] = 1000**2*theta*m_path_delta_geo[ti][3]
    # # extinct_mom_b_expon[ti] = 1 - (2*m_path_poisson[ti][0]**2) / (m_path_poisson[ti][2]) #+ scipy.stats.expon.pdf(x = 0, scale = 1/((2*m_path_poisson[ti][0])/(m_path_poisson[ti][2])))  
    # extinct_mom_b_power[ti] = 1 - (m_path_poisson[ti][0]**2)/(m_path_poisson[ti][2]) + scipy.stats.powerlaw.pdf(x=0, a = ((3*m_path_poisson[ti][2] - 2*m_path_poisson[ti][0])/(m_path_poisson[ti][2]-m_path_poisson[ti][0])))

# Geometric distribution reconstruction
# from scipy.stats import geom

# p = 1/shed_moments[0][100]
# x = np.arange(0,200)
# plt.plot(x, geom.pmf(x, p), 'bo', ms=8, label='geom pmf')
# plt.vlines(x, 0, geom.pmf(x, p), colors='b', lw=5, alpha=0.5)
# plt.show()

from datetime import datetime
date = datetime.today().strftime('%m-%d-%Y')

with open('extinction_mom_b_2layer_geometric_'+date+'_time500.npy', 'wb') as f:
    np.save(f, extinct_mom_b_geometric)
with open('extinction_mom_b_2layer_poisson_'+date+'_time500.npy', 'wb') as f:
    np.save(f, extinct_mom_b_poisson)
# with open('shed_first_moments_delta_02-26_time500_thetaalter.npy', 'wb') as f:
#     np.save(f, m_path_delta_geo[:,5])
# with open('basal_first_moment_geom_1-29_500.npy', 'wb') as f:
#     np.save(f, m_path_delta_geo[:,0])
# with open('para_first_moment_geom_1-29_500.npy', 'wb') as f:
#     np.save(f, m_path_delta_geo[:,1])

with open('para_second_moment_geom_'+date+'_time500.npy', 'wb') as f:
    np.save(f, m_path_delta_geo[:,3])

# # plt.plot(t_vec, cumu_extinct_poisson, label = 'Cumulative probability of extinction - Explicit basals (Poisson intital conditions)')
# plt.plot(t_vec, cumu_extinct_delta, label = 'Cumulative probability of extinction - Explicit basals (Delta intital conditions)')
# plt.plot(t_vec, extinct_mom_b_geometric, label = 'Cumulative probability of extinction - MOM basals (Geometric approx & Delta Initial conditions)')

# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# plt.title('2 Layer system')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()

# ## Parabasal plotting
# plt.plot(t_vec, cumu_extinct_delta_last, label = 'Cumulative probability of extinction - Explicit parabasals (Delta intital conditions)')
# plt.plot(t_vec, extinct_mom_p_geometric, marker = 'o', label = 'Cumulative probability of extinction - MOM parabasals (Geometric approx & Delta Initial conditions)')
# # plt.plot(t_vec, extinct_mom_p_poisson, label = "Cumulative probability of extinction - MOM parabasals (Poisson approx & poisson conditions)")
# # plt.plot(t_vec, cumu_extinct_poisson_last, label = "Cumulative probability of extinction - Explicit parabasals (Poisson approx & poisson conditions)")
# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Time')
# plt.title('2 Layer system')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()



# for time in range(0,200, 20):
#     plt.plot(np.arange(1,10000), scipy.stats.geom.pmf(k = np.arange(1,10000), p = 1/shed_moments[0][time]), label=fr"$t = {time}$")

# plt.legend()
# plt.ylabel('Probability')
# plt.xlabel('Virion count')
# plt.title('2 Layer system - Viral load')
# # plt.savefig("extinction_prob.pdf", format = 'pdf')
# plt.show()
# plt.close()










# with open('extinct_mom_b_2layer_main_geometric_1_8_500.npy', 'wb') as handle:    
#     np.save(handle, extinct_mom_b_geometric)
# with open('extinct_mom_b_2layer_main_delta_11_20.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_delta_geo)

# with open('extinct_mom_b_2layer_main_geometric_11_20.npy', 'wb') as handle:
#     np.save(handle, extinct_mom_b_geometric)
