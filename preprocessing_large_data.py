import random
import matplotlib.pyplot as plt
import numpy as np
import heapq
import math
import csv
import scipy



###################### Data Processing ##############################

# Input values
num_sims = 50000
tmax = 750
dt = 100

date = "04-26-2024"

with open('/Users/mcboudre/MOCS2/testCode/HPVCell/four_layer/data/basal_history_nonextinct_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_basal_history_nonextinct = np.load(f)

with open('/Users/mcboudre/MOCS2/testCode/HPVCell/four_layer/data/dead_history_all_'+date+'_time%d_sims%d.npy'%(tmax,num_sims), 'rb') as f:
    four_basal_history_nonextinct = np.load(f)
