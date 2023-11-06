"""Contains all user-defined event types.

All user-defined events are defined here, which
must be derived from the parent class EventBase.  

See also
---------
Module: base.events for documentation about the methods possible(), get_rate(), do_event() and get_invovle_other().

"""

import numpy as np
from base.events import EventBase

k_b = 8.617333262145e-5 # Boltzmann eV K-1

E_ads_A = -1.0 # A adsorption
E_ads_B = -1.0 # B adsorption 

E_act = 1.0 # Reaction barrier A+B

nu_des_A = 4e14  # Prefactor A desorption
nu_des_B2 = 1.0e10 # Prefactor B2 desorption
nu_react = 1.5e13  # Prefactor AB formation

class AAdsEvent(EventBase):
    """A adsorption event class.
    The event is A(g) + * -> A*.
    """

    def __init__(self, params):
        EventBase.__init__(self, params, name='AAds')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 2.0e2 * self.params['pA']
        return  R 

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 1 

    def get_involve_other(self):
        return False

class ADesEvent(EventBase):
    """A desorption event class.
    The event is A\* -> A + \*.
    """

    def __init__(self, params):
        EventBase.__init__(self, params,name='ADes')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = nu_des_A*np.exp(E_ads_A/(k_b*self.params['T'])) 
        return R 

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0

    def get_involve_other(self):
        return False

class B2AdsEvent(EventBase):
    """B2 adsorption event class.
    The event is 2B(g) + 2* -> 2B*.
    The event is possible if two neighbouring sites are empty.  
    The rate is set constant for comparison with mean-field model.  
    Performing the event adds two B to two sites.
    """

    def __init__(self, params):
        EventBase.__init__(self, params, name='B2Ads')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R =  4.0e3* self.params['pB2']
        return  R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 2
        system.sites[other_site].covered = 2

    def get_involve_other(self):
        return True

class B2DesEvent(EventBase):
    """A desorption event class.
    The event is 2B* -> B2(g) + 2*.
    The event is possible if two neighbouring sites are B-covered.  
    The rate comes from the forward rate and the
    equilibrium constant.  
    Performing the event removes two B from two sites.
    """

    def __init__(self, params):
        EventBase.__init__(self, params,name='B2Des')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 2 and system.sites[other_site].covered == 2:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = nu_des_B2*np.exp(2.*E_ads_B/(k_b*self.params['T']))
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

    def get_involve_other(self):
        return True

class ABreactEvent(EventBase):
    """A with B reaction event class.
    The event is A*+B* -> AB(g) + 2*.
    """

    def __init__(self, params):
        EventBase.__init__(self, params,name='ABreact')

    def possible(self, system, site, other_site):
        cov = [system.sites[site].covered, system.sites[other_site].covered] 
        if 1 in cov and 2 in cov:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = nu_react*np.exp(-E_act/(k_b*self.params['T'])) 
        return R  

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

    def get_involve_other(self):
        return True

class ADiffEvent(EventBase):
    """A diffusion event class.
    The event is A* + * -> * + A*.
    """

    def __init__(self, params):
        EventBase.__init__(self, params,name='ADiff')

    def possible(self, system, site, other_site):
        cov=[system.sites[site].covered, system.sites[other_site].covered] 
        if 0 in cov and 1 in cov:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1e8 
        return R 


    def do_event(self, system, site, other_site):
        old_cov_site = system.sites[site].covered
        old_cov_other_site = system.sites[other_site].covered
        system.sites[site].covered = old_cov_other_site 
        system.sites[other_site].covered = old_cov_site 

    def get_involve_other(self):
        return True 

class BDiffEvent(EventBase):
    """B diffusion event class.
    The event is B* + * -> * + B*.
    """

    def __init__(self, params):
        EventBase.__init__(self, params,name='BDiffEvent')

    def possible(self, system, site, other_site):
        if (system.sites[site].covered == 2 and system.sites[other_site].covered == 0) or (system.sites[site].covered == 0 and system.sites[other_site].covered == 2):
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1e8
        return R 


    def do_event(self, system, site, other_site):
        old_cov_site = system.sites[site].covered
        old_cov_other_site = system.sites[other_site].covered
        system.sites[site].covered = old_cov_other_site 
        system.sites[other_site].covered = old_cov_site 

    def get_involve_other(self):
        return True 

