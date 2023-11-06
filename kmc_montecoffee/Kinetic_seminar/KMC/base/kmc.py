"""Defines the NeighborKMCBase class.

The methods are used to perform kMC 
simulations with the first reaction method.

"""
from __future__ import print_function
from six.moves import configparser
import six

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser

import numpy as np
import random
random.seed()
#from bisect import bisect
import bisect
from sortedcontainers import SortedList

class NeighborKMCBase:
    """Main class for performing MonteCoffee simulations.
          
    Assigns a system to the simulation, stores parameters, 
    and reads in software configuration from the separate  
    file kMC_options.cfg.  
    Then it sets the time equal to zero and prepares to perform
    frm kinetic Monte Carlo simulations.

    Attributes:
    
    system: System
        The system instance to perform the simulation on.

    tend: float
        Simulation end-time, given in seconds.

    parameters: dict
        parameters used to calculate rate-constants and to dump to log files.
        Example: parameters = {'pCO':1E2,'T':700,'Note':'Test simulation'}

    t: float
        Current simulation time in seconds.

    *Attributes used to keep track of events (where and when they happen)*

    siteslist: list(int)
        The list of sites for each specific event.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    other_sitelist: list(int)
        The list of neighbor sites for each specific event.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    lastsel: int
        The int of the last selected site.

    lastother: int
        The int of the last selected neighbor site.

    rindex: list(list(list(int)))):
        The index of the specific events in lists like self.frm_times. For example to find the indices
        of site no i and event no j and neighbor number k to site i, call
        rindex[i][j][k].

    evs: numpy.ndarray(int):
        The event numbers for each specific event.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    rs: numpy.ndarray(float)
        Rate constants of specific events.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    wheres: list(list(int)):
        List of all the positions of the event-types in the lists with length
        len(self.events)*len(self.sites)*len(self.sites). To find all site-indices where event i
        happens, call wheres[i].
  
    involve_other: bool:
        -False if the event happens only on one specific site
        -True if the event modifies two or more sites         


    *Statistics counting attributes used to log and write output*

    SaveSteps: int
        The number of Monte Carlo steps between saving the .txt files.

    LogSteps: int
        The number of Monte Carlo steps between logging steps.

    tinfinity: float
        What time to put impossible events to.

    Nspecies: int
        How many different types of species are in the simulation. Used to
        print and log.
   
    nninter: int
        How deep is the nearest-neighbor interaction (depth of effect of event on neighbor properties)

    verbose: bool
        If True, the code prints verbose information.

    save_coverages: bool
        If True, coverages are saved to coverages.txt and the site, othersite and event evolution in detail_site_event_evol.hdf5. This can result in
        large files.

    write_atoms: bool
        If True, the surface atoms are written with the step number in the filename. It has to be adjusted for adsorption species individually. 

    times: list(float)
        List of times for each logged monte carlo steps in self.MCstep

    MCstep: list(int)
        List of Monte Carlo step numbers logged.

    Nsites: int
        The number of sites in self.system

    Nstypes: int
        The number of distinct site-types.

    covered: list(list(int))
        A list of site-occupations, of each site for each logged step.
        To find the site-occupation at step no self.MCstep[i] and site j, call
        covered[i][j].

    system_evolution: list(list())
        List which contains a list of the fired event with at site, other site and time

    used_ijk: list(tuples(site,event,othersite))
        List of tuples representing the unique neighbor-event pairs avoiding double counting. 


    """

    def __init__(self, system, tend, parameters={}, tstart = 0):

        self.system = system
        self.tend = tend
        self.parameters = parameters

        self.t =  tstart  #Initialize the time

        self.update_all_rates = False # Shall we update all the events or just a part?!
        # Load software configuration
        self.load_options()

        if self.verbose:
            print('-' * 50, '\n', 'MonteCoffee Simulation Initialized', '\n', '-' * 50, '\n')
            print('kMC simulation loading ...')

        # Variables connected to after analysis.
        self.times = []
        self.covered = []  # Site-occupations
        self.evs_exec = np.zeros(len(self.events))
        self.system_evolution = [[] for i in range(4)]

        self.Nstypes = list(set([s.stype for s in self.system.sites]))

        #self.used_ijk = [] # Avoid double counting events of site 1+2 vs site 2+1.

        # FRM method variables
        # --------------------------------------------------     
        self.frm_times = []  # Times of occourences
        self.frm_arg = None  # args that sort frm times
        if self.verbose:
            print('Initializing First Reaction method lists ...')

        self.frm_init()

    def load_options(self):
        """Loads all options set in kMC_options.cfg.
        
        Instantiates a configuration parser, and loads in all
        options from *kMC_options.cfg*.

        """
        config = configparser.RawConfigParser()
        config.read('kMC_options.cfg')
        self.SaveSteps = config.getint('Parameters', 'SaveSteps', fallback = 1000)  # How often to save txt files
        self.LogSteps = config.getint('Parameters', 'LogSteps', fallback = 1)  # How often to log steps
        self.tinfinity = config.getfloat('Parameters', 'tinfinity', fallback = 1e18)  # What is considered infinite time
        self.Nspecies = config.getint('Parameters', 'Nspecies', fallback = 1)  # Number of different species in simulation.
        self.nninter = config.getint('Parameters', 'nninteractions', fallback = 1)  # depth of NN interactions
        self.verbose = config.getboolean('Options', 'Verbose', fallback = True)  # Print verbose information?
        self.write_atoms = config.getboolean('Options', 'Write_atoms', fallback = False)  # Write out the atom object
        self.save_coverages = config.getboolean('Options', 'SaveCovs', fallback = False)  # Save coverages?

    def frm_init(self):
        """Prepare to perform FRM simulation.
            
        Initializes empty rate and event lists to 
        bookkeep the FRM algorithm. The initial times
        of occurrence for each event is also calculated
        and stored.  

        """
        self.lastsel = 0
        self.lastother = None
        self.possible_ijk = set()    # possible tubles of all possible events in same order as self.possible_times
        self.possible_times = SortedList()    # Times of all possible events
        l_used_ijk = set()         # Employed site_event_othersite paires (unique) 

        for j, e in enumerate(self.events):
            for i, s in enumerate(self.system.sites):
                NNcur = self.system.neighbors[i]
                if not e.get_involve_other():
                    if e.possible(self.system, i, NNcur[0]): 
                        rcur = e.get_rate(self.system, i, NNcur[0])
                        u = random.uniform(0,1)
                        self.possible_times.add((self.t - np.log(u) / rcur,(i,j,NNcur[0])))
                        self.possible_ijk.add((i,j,NNcur[0]))
                    l_used_ijk.add((i,j,NNcur[0])) 
                else:
                    for k, other_site in enumerate(NNcur):
                        used_ijk = True if (i,j,other_site) in l_used_ijk or (other_site,j,i) in l_used_ijk else False # Check if the pair is already in list 
                        if e.possible(self.system, i, other_site) and not used_ijk: 
                            rcur = e.get_rate(self.system, i, other_site)
                            u = random.uniform(0,1)
                            self.possible_times.add((self.t - np.log(u) / rcur,(i,j,other_site)))
                            self.possible_ijk.add((i,j,other_site))
    
                        if (i,j,other_site) not in l_used_ijk and (other_site,j,i) not in l_used_ijk: # Create unique pair-event tuple list  
                            l_used_ijk.add((i,j,other_site)) 

#        self.possible_times, self.possible_ijk = map(list, zip(*sorted(zip(self.possible_times,self.possible_ijk))))
#        self.possible_times = SortedList(self.possible_times)
#        self.used_ijk = set(self.used_ijk)
#        self.frm_arg = np.argmin(self.possible_times)

    def find_index(self, check_tuble):
        for ind, element in enumerate(self.possible_times):
            if check_tuble == element[1]:
                return ind               


    def frm_update(self):

        search = self.system.find_nn_recurse(self, [self.lastsel,self.lastother])  # returns a sorted list
        updated_set = set()
        for j, e in enumerate(self.events):
            if not e.get_involve_other():
                for i in search:
                    NNcur = self.system.neighbors[i][0]
                    poss_now = e.possible(self.system, i, NNcur)
                    in_poss_list = True if (i,j,NNcur) in self.possible_ijk else False
                    if poss_now and not in_poss_list:
                        rcur = e.get_rate(self.system, i, NNcur)
                        u = random.uniform(0,1)
                        new_time = self.t - np.log(u) / rcur
                        self.possible_times.add((new_time,(i,j,NNcur))) 
                        self.possible_ijk.add((i,j,NNcur))
                    if not poss_now and in_poss_list:
                        poss_index = self.find_index((i,j,NNcur))
                        self.possible_ijk.remove((i,j,NNcur))
                        del self.possible_times[poss_index]

            else:
                for i in search:
                    for other in self.system.neighbors[i]:
                        i_ass = min([i,other])
                        other_ass = max([i,other])
                        updated = True if (i_ass,j,other_ass) in updated_set else False 
                        if not updated:
                            poss_now = e.possible(self.system, i_ass, other_ass)
                            in_poss_list = True if (i_ass,j,other_ass) in self.possible_ijk else False
                            if poss_now and not in_poss_list:
                                rcur = e.get_rate(self.system, i_ass, other_ass)
                                u = random.uniform(0,1)
                                new_time = self.t - np.log(u) / rcur
                                self.possible_times.add((new_time,(i_ass,j,other_ass))) 
                                self.possible_ijk.add((i_ass,j,other_ass))
        
                            if not poss_now and in_poss_list:
                                poss_index = self.find_index((i_ass,j,other_ass))
                                self.possible_ijk.remove((i_ass,j,other_ass))
                                del self.possible_times[poss_index]

                            updated_set.add((i_ass,j,other_ass))

        # Update the first reaction part
#        self.frm_arg = np.argmin(self.possible_times)

    def frm_step(self):

        # Choose the first reaction if possible
        site = self.possible_times[0][1][0]  # The site to do event.
        othersite = self.possible_times[0][1][2]
        self.lastsel = int(site)
        self.lastother = int(othersite)
        evtype = self.possible_times[0][1][1]
        if self.events[evtype].possible(self.system, site, othersite):
            # Event is possible, change state
            self.events[evtype].do_event(self.system, site, othersite)
            self.t = self.possible_times[0][0]  # Update time

            #print ('Event', evtype, self.possible_times[0])
            #Data logging
            self.evs_exec[evtype] += 1 #Total executed events
            self.system_evolution[0].append(int(self.system.sites[site].ind))
            self.system_evolution[1].append(int(self.system.sites[othersite].ind))
            self.system_evolution[2].append(int(evtype))
            self.system_evolution[3].append(float(self.t))
            
            #Delete the executed event from possible list
            self.possible_ijk.remove((site,evtype,othersite))
            del self.possible_times[0]

        else:
#            del self.possible_times[self.frm_arg]
#            del self.possible_ijk[self.frm_arg]
            # New first reaction must be determined
            print ('Event attempted', evtype, self.possible_times[0])
            print (self.system.sites[site].covered, self.system.sites[othersite].covered)
            raise Warning("Impossible event were next in que and was attempted")

        self.frm_update()

    def update_T(self, T_init):
        heat_rate = self.parameters['heating_rate']
        upd_T = T_init + self.t*heat_rate
        self.parameters['T'] = upd_T
        self.update_all_rates = True

 
    def load_events(self):
        """Loads events (abstract method).

        This method must be overridden by the child class in user_kmc.NeighborKMC.

        Raises
        ---------
        NotImplementedError:
            If called.
          
        """
        raise NotImplementedError('''User needs to define load_events
                                 method in derived NeighborKMC class''')

    def run_kmc(self):
        """Runs the kMC simulation (abstract method)

        This method must be overridden by the child class in user_kmc.NeighborKMC.

        Raises
        ---------
        NotImplementedError:
            If called.

        """
        raise NotImplementedError('''User needs to define run_kmc method 
                                          in derived NeighborKMC class''')
