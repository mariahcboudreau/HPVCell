'''
Defines the stochastic cell model for HPV (Markov chain)
'''

import numpy as np
import random
import matplotlib.pyplot as plt
from . import cell_mass as cellMass


'''
TO DO:
Allow for random viral load to be floating around and be in the system
Normal cell division should happen, then allow there to be an infected section 
'''

'''
Event driven functions

Code that defines the functions that allow events to be drawn at varying time points over continuous time. 
'''
class Event_Driven:


    def draw_event_basal_normal(self):
        '''
        Draws the type of event that could occur for each basal cell (division, or infection)

        Returns:
            (int): event class value
        '''

        basal_normal_bb_rate = 0  # draws the event class possibility
        basal_normal_pp_rate = 0
        basal_normal_bp_rate = 0
        infect_rate = 0


        # Construct the vectors
        basal_inds_normal = []

        for i in basal_inds_normal:
            basal_normal_bp_rate += cellPar.get_division_rate(i) #TODO access the split rates in parameters
            basal_normal_bb_rate += cellPar.get_division_rate(i)
            basal_normal_pp_rate += cellPar.get_division_rate(i)
            infect_rate += cellPar.get_infect_rate(i)



        basal_normal_bp_start = 0
        basal_normal_bp_end = basal_normal_bp_rate
        basal_normal_bb_end = basal_normal_bp_end + basal_normal_bb_rate
        basal_normal_pp_end = basal_normal_bb_end + basal_normal_pp_rate
        infect_end = basal_normal_pp_end + infect_rate



        random_draw = random.uniform(basal_normal_bp_start, infect_end)

        if random_draw < basal_normal_bp_end:
            return 4                       # asymmetric normal split (BP) from basal
        elif (random_draw >= basal_normal_bp_end) & (random_draw < basal_normal_bb_end):
            return 3                        # symmetric normal split (BB) from basal
        elif (random_draw >= basal_normal_bb_end) & (random_draw < basal_normal_pp_end):
            return 2                        # symmetric normal split (PP) from basal
        elif (random_draw >= basal_normal_pp_end) & (random_draw < infect_end):
            return 1                        # infection event


    def draw_event_basal_infect(self):
        '''
        Draws the type of event that could occur for each basal cell (division, or transformation)

        Returns:
            (int): event class value
        '''


        basal_infect_bb_rate = 0
        basal_infect_pp_rate = 0
        basal_infect_bp_rate = 0
        transform_rate = 0

        # Construct the vectors TODO filter here
        basal_inds_infect = []


        for i in basal_inds_infect:
            basal_infect_bb_rate += cellPar.get_division_rate(i)
            basal_infect_pp_rate += cellPar.get_division_rate(i)
            basal_infect_bp_rate += cellPar.get_division_rate(i)
            transform_rate += cellPar.get_transform_rate(i)

        basal_infect_bp_start = 0
        basal_infect_bb_end = basal_infect_bp_start + basal_infect_bb_rate
        basal_infect_bp_end = basal_infect_bb_end + basal_infect_bp_rate
        basal_infect_pp_end = basal_infect_bp_end + basal_infect_pp_rate
        transform_end = basal_infect_pp_end + transform_rate


        random_draw = random.uniform(basal_infect_bp_start, transform_end)

        if random_draw < basal_infect_bb_end:
            return 4                        # symmetric infected split (BB) from basal
        elif (random_draw >= basal_infect_bb_end) & (random_draw < basal_infect_bp_end):
            return 3                        # asymmetric infected split (BP) from basal
        elif (random_draw >= basal_infect_bp_end) & (random_draw < basal_infect_pp_end):
            return 2                        # symmetric infected split (PP) from basal
        elif(random_draw >= basal_infect_pp_end) & (random_draw < transform_end):
            return 1                        # transformation event


    def draw_event_parabasal_normal(self):
        '''
        Draws the type of event that could occur for each parabasal cell (division, or differentiation)


        Updates the Cell Mass Array

        '''


        pbasal_normal_pp_rate = 0
        diff_normal_rate = 0

        # Construct the vectors
        pbasal_inds_normal = [] #TODO alter these for filtering

        for i in pbasal_inds_normal:
            pbasal_normal_pp_rate += cellPar.get_division_rate(i) #TODO access the split rates in parameters
            diff_normal_rate += cellPar.get_diff_rate(i)


        pbasal_normal_pp_start = 0
        pbasal_normal_pp_end = pbasal_normal_pp_rate
        diff_normal_end = pbasal_normal_pp_end + diff_normal_rate


        random_draw = random.uniform(pbasal_normal_pp_start, diff_normal_end)

        if random_draw < pbasal_normal_pp_start:
            # Assign 5 to the attribute for that cell
            return 2                      # symmetric normal split (PP) from parabasal
        elif (random_draw >= pbasal_normal_pp_start) & (random_draw < diff_normal_end):
            return 1




    def draw_event_parabasal_infect(self):
        '''
        Draws the type of event that could occur for each parabasal cell (division, differentiation or transformation)


        Updates the Cell Mass Array

        '''


        pbasal_normal_pp_rate = 0
        pbasal_infect_pp_rate = 0
        transform_rate = 0
        diff_normal_rate = 0
        diff_infect_rate = 0

        # Construct the vectors
        pbasal_inds_normal = [] #TODO alter these for filtering
        pbasal_inds_infect = []

        for i in pbasal_inds_normal:
            pbasal_normal_pp_rate += cellPar.get_division_rate(i) #TODO access the split rates in parameters
            diff_normal_rate += cellPar.get_diff_rate(i)

        for i in pbasal_inds_infect:
            pbasal_infect_pp_rate += cellPar.get_division_rate(i)
            diff_infect_rate += cellPar.get_diff_rate(i)
            transform_rate += cellPar.get_transform_rate(i)

        pbasal_normal_pp_start = 0
        pbasal_normal_pp_end = pbasal_normal_pp_rate
        pbasal_infect_pp_end = pbasal_normal_pp_end + pbasal_infect_pp_rate
        transform_end = pbasal_infect_pp_end + transform_rate
        diff_normal_end = transform_end + diff_normal_rate
        diff_infect_end = diff_normal_end + diff_infect_rate

        all_inds = pbasal_inds_normal + pbasal_inds_infect
        for i in all_inds:
            random_draw = random.uniform(pbasal_normal_pp_start, diff_infect_end)

            if random_draw < pbasal_normal_pp_start:
                # Assign 5 to the attribute for that cell
                return 5                        # symmetric normal split (PP) from parabasal
            elif (random_draw >= pbasal_normal_pp_start) & (random_draw < pbasal_infect_pp_end):
                return 4                        # symmetric infected split (PP) from parabasal
            elif (random_draw >= pbasal_infect_pp_end) & (random_draw < transform_end):
                return 3                        # transformation event
            elif (random_draw >= transform_end) & (random_draw < diff_normal_end):
                return 2                        # differentiation event on a normal parabasal
            elif (random_draw >= diff_normal_end) & (random_draw < diff_infect_end):
                return 1                        # differentiation event on an infected parabasal


def create_cell_events():
    '''


    Returns:
        parabasal event
    '''


def draw_event_class(self, V_B, V_P_non_diff, V_infected):
        '''
        Draws the type of event that will could happen, either cell splitting (basal cell to two new basal cells, two new parabasal cells, or a parabasal cell and a basal cell,
        or parabasal cell to two new parabasal cells,
        or parabasal cells differentiating,
        or basal cell becoming infected )

        Args:
            V_B (list/basal cells): vector of basal cells
            V_P_non_diff (list/parabasal cells): vector of parabasal cells (have not differentiated yet)
            V_infected (array/infected basal cells and infected parabasal cells): vector of infected cells

        Returns:
            (int): event class value

        '''
        basal_split_bb_rate = 0                     # draws the event class possibility
        basal_split_pp_rate = 0
        basal_split_bp_rate = 0
        pbasal_split_pp_rate = 0
        infect_rate = 0
        diff_rate = 0
        transform_rate = 0

        for cell in V_B: # TODO: Make sure that this is correct and do not need to have different event rates
            basal_split_bp_rate += cell.event_rate
            basal_split_bb_rate += cell.event_rate
            basal_split_pp_rate += cell.event_rate
            infect_rate += cell.event_rate

        for cell in V_P_non_diff:
            pbasal_split_pp_rate += cell.event_rate
            diff_rate += cell.event_rate

        for cell in V_infected:
            transform_rate += cell.event_rate


        basal_bp_start = 0
        basal_bp_end = basal_split_bp_rate
        basal_bb_end = basal_bp_end + basal_split_bb_rate
        basal_pp_end = basal_bb_end + basal_split_pp_rate
        pbasal_pp_end = basal_pp_end + pbasal_split_pp_rate
        infect_end = pbasal_pp_end + infect_rate
        diff_end = infect_end + diff_rate
        trans_end = diff_end + transform_rate


        random_draw = random.uniform(basal_bp_start, trans_end)

        if random_draw < basal_bp_end:
            return 7                                                        # asymmetric split (BP) from basal
        elif (random_draw >= basal_bp_end) & (random_draw < basal_bb_end):
            return 6                                                        # symmetric split (BB) from basal
        elif (random_draw >= basal_bb_end) & (random_draw < basal_pp_end):
            return 5                                                        # symmetric split (PP) from basal
        elif (random_draw >= basal_pp_end) & (random_draw < pbasal_pp_end):
            return 4                                                        # symmetric split (PP) from parabasal
        elif (random_draw >= pbasal_pp_end) & (random_draw < infect_end):
            return 3                                                        # infection event
        elif (random_draw >= infect_end) & (random_draw < diff_end):
            return 2                                                        # differentiation event
        elif (random_draw >= diff_end) & (random_draw < transform_rate):
            return 1                                                        # transformation event


    def draw_event(max_rate, event_list):
        '''
        Draws the type of event that will could happen, either cell splitting (basal cell to two new basal cells, two new parabasal cells, or a parabasal cell and a basal cell,
        or parabasal cell to two new parabasal cells,
        or parabasal cells differentiating,
        or basal cell becoming infected,
        or infected cell becoming transformed)

        Args:
            max_rate (float): largest rate in the vector
            event_list (list/per event): vector of cells that could have an event accepted to happen

        Returns:
            cell (Cell): cell that has an event happening to it
        '''
        accepted = False
        random_event = None
        while not accepted:
            random_event = random.choice(event_list)
            accept_rate = random_event.event_rate / max_rate
            random_draw = random.uniform(0,1)
            if random_draw < accept_rate:
                accepted = True

        return random_event

'''
CLASS DECLARATIONS
'''
class MarkovSim:                                    # Simulation to let the event-driven stochastic process occur
    '''
    Simulation to let the event-driven stochastic cell cycle occur
    Args:
        time (int): total time the process will run
        b_bb (double): rate informing symmetric split of basal cell to two new basal cells
        b_pb (double): rate informing asymmetric split of basal cell to a basal cell and parabasal cell
        b_pp (double): rate informing symmetric split of basal cell to two new parabasal cells
        p_pp (double): rate informing symmetric split of parabasal cell to two new parabasal cells
        diff (double): rate informing differentiation of parabasal cells
        infect_cell (double): rate informing the infection of basal cells through microabrasions

    '''

    def __init__(self, time, b_bb, b_pb, b_pp, p_pp, diff, infect_cell, transform):
        self.indices = list(range(1, 1000))         # number of cells in the system
        self.start_num_cells = list(range(1, 900))  # number of cells in the system at the start
        self.time = time                            # total time of simulation
        self.current_time = 0                       # current time in simulation; starts at 0
        self.b_bb = b_bb                            # rate informing basal cell symmetric split to basal cells
        self.b_pb = b_pb                            # rate informing basal cell asymmetric split
        self.b_pp = b_pp                            # rate informing basal cell symmetric split to parabasal cells
        self.p_pp = p_pp                            # rate informing parabasal cell symmetric split
        self.diff = diff                            # rate informing differentiation
        self.transform = transform                  # rate informing transformation of infected cell
        self.infect_cell = infect_cell              # rate informing infection
        self.V_B = []                               # vec of basal cells (objects)
        self.V_B_ind = []                           # vec of basal cell indices
        self.V_P_non_diff = []                      # vec of parabasal cells that are not differentiating
        self.V_P_non_diff_ind = []                  # vec of indices of above
        self.V_P_non_diff_infect = []               # vec of parabasal cells that are not differentiating and infected
        self.V_P_non_diff_infect_ind = []           # vec of indices of above
        self.V_P_diff = []                          # vec of parabasal cells that have potential to divide
        self.V_P_diff_ind = []                      # indices of potential dividing parabasal cells
        self.vl_time = np.zeros(time)               # measured in weeks?  TODO: determine a time interval
        self.V_infected = []                    # vector for infected cells
        self.V_infected_indices = []            # indices of infected cells
        self.num_infected_t = np.zeros(self.time)   # counter for the cells infected at each time point TODO: change this as cells are infected
        self.deaths = []                            # times at which a differentiated cell will shed and die
        self.shed_amount_t =  np.zeros(self.time)   # counter for the amount of viral load shed at each time point

    def initialize(self): # TODO: initialize 6% of cells to be basal, the rest parabasal and 75% of parabasal are diff
        base_cells = []
        for i in range(self.start_num_cells):
            base_cells.append(Cell(i, 1,1))
            if i < (.06 * len(self.indices)):
                self.V_B.append(Basal_Cell(base_cells[i].index, base_cells[i].split_rate, base_cells[i].death_rate, base_cells[i]))
            else:
                if i > (0.705 * len(self.indices)):
                    self.V_P_diff.append(Parabasal_Cell(base_cells[i].index, base_cells[i].split_rate, base_cells[i].death_rate, base_cells[i]))
                else:
                    self.V_P_non_diff.append(Parabasal_Cell(base_cells[i].index, base_cells[i].split_rate, base_cells[i].death_rate,
                                       base_cells[i]))


    def run_sim(self):
        '''
        Runs one iteration of cell cycle and possible infection events
        '''
        self.initialize()
        while self.current_time < self.time:

            tau = 1
            event_class = draw_event_class(self.V_B, self.V_P_non_diff)
            index_1 = np.random.choice(self.indices)
            index_2 = np.random.choice(self.indices)
            while index_2 == index_1:
                index_2 = np.random.choice(self.indices)

            if event_class == 1:
                # Transforming an infected cell
                trans_event = draw_event(np.max(self.transform), self.V_infected)
                trans_event.transform()

                #Bookkeeping of transformation and how that affects cell death/split rate

            if event_class == 2:
                # Differentiating a parabasal cell
                diff_event = draw_event(np.max(self.Lambda), self.V_P_non_diff)
                diff_event.differentiate()

                # Bookkeeping for differentiating parabasal cells
                self.V_P_diff.append(diff_event)
                self.V_P_diff_ind.append(diff_event.index)
                self.V_P_non_diff.remove(diff_event)
                self.V_P_non_diff_ind.remove(diff_event.index)
                self.deaths.append(diff_event.death_time)           # appends a time when the cells die



            if event_class == 3:
                # Infect a basal cell
                infect_event = draw_event(np.max(self.Lambda), self.V_B)
                infection = Infected_Basal_Cell('16', 500, infect_event)
                self.infected_cells.append(infection)
                self.infected_cells_indices.append(infection.cell.index)

            if event_class == 4:
                # Make two new parabasal cells from parabasal cell
                pbasal_pp_event = draw_event(np.max(self.Lambda), self.V_P_non_diff_infect)
                new_p1, current = pbasal_pp_event.split(index_1)

                #Bookkeeping and updating the lists of cells
                self.V_P_non_diff.append(new_p1)
                self.V_P_non_diff_ind.append(new_p1.index)
                self.indices.remove(new_p1.index)

                #Bookkeeping the infection
                if isinstance(new_p1, Infected_Parabasal_Cell):
                    self.infected_cells.append(new_p1)
                    self.infected_cells_indices.append(new_p1.index)
                    if new_p1.state == 0:
                        self.V_P_non_diff_infect.append(new_p1)
                        self.V_P_non_diff_infect_ind.append(new_p1.index)


            if event_class == 5:
                # Make two new parabasal cells from basal cell
                basal_pp_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p1, new_p2 = basal_pp_event.split(4, index_1, index_2, self.current_time)

                #Bookkeeping and updating the lists of cells
                self.V_P_non_diff.extend((new_p1, new_p2))
                self.V_P_non_diff_ind.extend((new_p1.index, new_p2.index))
                self.indices.remove(new_p1.index)
                self.indices.remove(new_p2.index)
                self.V_B.remove(basal_pp_event)
                self.indices.append(basal_pp_event.index)

                # Bookkeeping the infection
                if isinstance(new_p1, Infected_Parabasal_Cell):
                    self.infected_cells.append(new_p1)
                    self.infected_cells_indices.append(new_p1.index)
                    if new_p1.state == 0:
                        self.V_P_non_diff_infect.extend((new_p1, new_p2))
                        self.V_P_non_diff_infect_ind.extend((new_p1.index,new_p2.index))

            if event_class == 6:
                # Make two new basal cells from basal cell
                basal_bb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_b1, current = basal_bb_event.split(5, index_1, basal_bb_event.index, self.current_time)

                # Bookkeeping and updating the lists of cells
                self.V_B.append(new_b1)
                self.V_B_ind.append(new_b1.index)
                self.indices.remove(new_b1.index)

                # Bookkeeping the infection
                if isinstance(new_b1, Infected_Basal_Cell):
                    self.infected_cells.append(new_b1)
                    self.infected_cells_indices.append(new_b1.index)

            if event_class == 7:
                # Make a basal and a parabasal cell from basal cell
                basal_pb_event = draw_event(np.max(self.Lambda), self.V_B)
                new_p, current = basal_pb_event.split(6, index_1, basal_pb_event.index, self.current_time)

                # Bookkeeping and updating the lists of  cells
                self.V_P_non_diff.append(new_p)
                self.V_P_non_diff_ind.append(new_p.index)
                self.indices.remove(new_p.index)

                # Bookkeeping the infection
                if isinstance(new_p, Infected_Parabasal_Cell):
                    self.infected_cells.append(new_p)
                    self.infected_cells_indices.append(new_p)
                    if new_p.state == 0:
                        self.V_P_non_diff_infect.append(new_p)
                        self.V_P_non_diff_infect_ind.append(new_p.index)



            for i in range(tau):
                if (self.current_time == self.deaths[0]):
                    self.deaths.pop(0)
                    if isinstance(self.V_P_diff[0], Infected_Parabasal_Cell):
                        self.shed_amount_t[self.current_time] += self.V_P_diff[0].shed()
                        self.shed_amount_t[self.current_time + 1] += self.V_P_diff[0].shed()/2
                        self.shed_amount_t[self.current_time + 2] += self.V_P_diff[0].shed()/4

                    self.V_P_diff[0].die()

                self.current_time += 1







#
# class Cell:
#     '''
#     Base class to give general attributes of a cell
#
#     Args:
#         index (int): index of the cell
#         split_rate (float): rate of division for the cell
#     '''
#
#
#     def __init__(self, index, split_rate, death_rate):
#         self.index = index
#         self.split_rate = split_rate
#         self.death_rate = death_rate
#
#     def update_split(self, new_rate): # E6/E7 expressions ?
#         self.split_rate = new_rate
#
#
#
# class Basal_Cell(Cell):
#     '''
#         Basal cell that can be infected and does most of the cell division
#
#         Args:
#             index (int): index of the cell
#             split_rate (float): rate of division for the cell
#             cell (Cell): the original cell that created this instance
#     '''
#     def __init__(self, index, split_rate, cell = None):
#         super().__init__(index, split_rate)
#         self.cell = cell
#
#     def split(self, event_type, index_1, index_2, time):
#         # Make the extra cell, default parabasal cell which is infected
#         if event_type == 5: # two parabasal cells
#             split_1 = Parabasal_Cell(index_1, self.split_rate, self, time)
#             split_2 = Parabasal_Cell(index_2, self.split_rate, self, time)
#
#             if isinstance(self, Infected_Basal_Cell):
#                 split_1 = Infected_Parabasal_Cell(self.genotype, self.vl, self, time)
#                 split_2 = Infected_Parabasal_Cell(self.genotype, self.vl, self)
#
#         if event_type == 6: # two new basal cells
#             split_1 = Basal_Cell(index_1, self.split_rate, self, time)
#             split_2 = self
#
#             if isinstance(self, Infected_Basal_Cell):
#                 split_1 = Infected_Basal_Cell(self.genotype, self.vl, self)
#
#         if event_type == 7: # one basal, one parabasal
#             split_1 = Parabasal_Cell(index_1, self.split_rate, self, time)
#             split_2 = self
#
#             if isinstance(self, Infected_Basal_Cell):
#                 split_1 = Infected_Parabasal_Cell(self.genotype, self.vl, self)
#
#         return split_1, split_2
#
#
#
# class Parabasal_Cell(Cell):
#     '''
#         Parabasal cell that can divide
#
#            Args:
#             index (int): index of the cell
#             split_rate (float): rate of division for the cell
#             cell (Cell): the original cell that created this instance
#             death_time (float): time the cells dies and potentially sheds
#             dead (bool): if cell is dead or not
#     '''
#
#
#     def __init__(self, index, split_rate, cell, time):
#         super().__init__(index, split_rate)
#         self.state = 0
#         self.cell = cell
#         self.death_time = time+3
#         self.dead = False
#
#     def differentiate(self):
#         self.state = 1
#
#     def is_dead(self):
#         return self.dead
#
#     def die(self):
#         self.dead = True
#
# class Infected_Basal_Cell(Basal_Cell):
#     '''
#         Infected basal cell
#
#         Args:
#             index (int): index of the cell
#             split_rate (float): rate of division for the cell
#             cell (Cell): the original cell that created this instance
#             death_time (float): time the cells dies and potentially sheds
#             dead (bool): if cell is dead or not
#         '''
#     # cell: the cell that has just become infected
#     # genotype: the genotype of the virus in the system
#     # vl: viral load of the cell
#     # transformed: transformed cell or not
#
#     def __init__(self, genotype, vl, cell):
#         super().__init__(cell.index, cell.split_rate, cell.death_rate, cell)
#         self.genotype = genotype
#         self.vl = vl
#         self.transformed = False
#
#     def update_vl(self):
#         self.vl = 1000
#
#     def shed(self):
#         return self.vl
#
#     def transform(self):
#         self.transformed = True
#         #change the death rate and cell splitting
#
#
#
# class Infected_Parabasal_Cell(Parabasal_Cell):
#     # cells: associated cells that the virus is inhabiting
#     # cell: the cell that has just become infected
#     # genotype: the genotype of the virus in the system
#     # vl: viral load of the cell
#     # shed_time: time that the cell sheds it viral load
#     # transformed: transformed cell or not
#
#     def __init__(self, genotype, vl, cell):
#         super().__init__(cell.index, cell.split_rate, cell.death_rate, cell.state, cell)
#         self.genotype = genotype
#         self.vl = vl
#         self.transformed = False
#
#     def split(self, index_1):
#         # Make the extra infected cell
#         split_1 = Parabasal_Cell(index_1, 0, self)
#         split_2 = self
#
#         if isinstance(self, Infected_Parabasal_Cell):
#             split_1 = Infected_Parabasal_Cell(self.genotype, self.vl, self)
#
#         return split_1, split_2
#
#     def update_vl(self):
#         self.vl = 1000
#
#     def shed(self):
#         return self.vl
#
#     def transform(self):
#         self.transformed = True
#         # change the death rate and cell splitting
#
