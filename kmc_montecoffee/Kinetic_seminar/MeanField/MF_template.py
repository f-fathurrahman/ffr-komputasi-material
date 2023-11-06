'''Import needed modules'''
import numpy as np
from scipy import integrate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

'''Reaction conditions & Physical constants'''
T = 450 # K
P_Ag =  # Pa
P_B2g =  # Pa
k_b = 8.617333262145e-5 # Boltzmann eV K-1
theta0 =  # set initial coverages to zero
'''Adsorption energies & Barrier'''

# A adsorption
# B2 adsorption

# A+B surface reaction barrier

'''Pre-exponential factors'''
# A desorption
# B2 desorption

# surface reaction

def get_rate_constants(T):
    '''This function returns forward and backwards rate constants at T'''

    kbT = k_b*T

    '''Adsorption rate constants'''
    # Adsorption of A
    # Adsorption of B2

    '''Desorption rate constants'''
    # Desorption of A
    # Desorption of B2

    '''Surface reaction rate constant'''
    # surface reaction

    return k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf

def get_rates(theta,k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf):
    '''This function returns reaction rates calculated from
    coverages and rate constants'''

    # coverage of A
    # coverage of B

    #site balance

    rate = [0]*3

    rate[0] = # Adsorption rate of A
    rate[1] = # Adsorption rate of B2
    rate[2] =  # Surface reaction rate

    return rate

def get_system_of_ODEs(t,theta,k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf):
    '''Returns differential equations'''

    rate = get_rates(theta,k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf)

    dthetadt = [0]*2
    dthetadt[0] =# dA/dt
    dthetadt[1] =# dB/dt

    return dthetadt #this HAS to be a list

def ode_solver(theta0,k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf):
    '''Solves set of ODEs for a given kf,kr, and initial thetas.
    I suggest you do not change this part
    unless you know what you're doing'''
    t_start = 0 #starting time
    t_final = 100000 #ending time
    steps = 10000 #sampling steps
    delta_t = (t_final - t_start) / steps
    t_eval = np.arange(t_start, t_final, delta_t)# Times solution is stored
    tspan = (t_start,t_final) # Interval of integration in seconds

    solution = solve_ivp(get_system_of_ODEs, # fun(t, y)
                    tspan,
                    theta0, # Initial state i.e. coverages
                    method = 'BDF', #backward differentiation formula based method
                    t_eval= t_eval,
                    #rtol=1e-10, # relative tolerance
                    #atol=1e-11, # absolute tolerance
                    args=(k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf)) #get_system_of_ODEs also takes arguments kf,kr
    print('-------------------------------------')
    print(solution.message)
    print('-------------------------------------')
    return solution.t, solution.y.T

def solve_microkinetic_model(T,theta0):

        k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf = get_rate_constants(T)

        '''solve mk model for calculated rate constants
            return time and coverages'''
        timepoints, coverages = ode_solver(theta0,k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf)

        return timepoints, coverages

def get_results(theta_f,T):
    '''Returns reaction rates for a given theta'''
    k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf = get_rate_constants(T)
    finalrates = get_rates(theta_f,k_ads_A,k_ads_B2,k_des_A,k_des_B2,k_surf)
    return finalrates

'''Solve ODEs for given temperature(s)'''

times,c = solve_microkinetic_model(T,theta0)
rate = get_results(c[-1],T)

'''Plot results'''
