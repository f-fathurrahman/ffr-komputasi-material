import matplotlib.pyplot as plt
import numpy as np

def add_boundaries(input_array):    
    '''adding extra voxels left and right to build the hot plates.
    adding extra voxels at front, back, top and bottom with kappa = 0 to model adiabatic conditions.'''
    length = input_array.shape[0]
    width = input_array.shape[1]
    height = input_array.shape[2]   
    # generate voxels with kappa = inf to model the hotplates for the left and right side
    new_voxels_lr = np.ones((1, width, height)) * float('inf') 
    # generate voxels with kappa = 0 to model the hotplates for top, bottom, front and back
    new_voxels_bf = np.zeros((length+2, 1, height)) 
    new_voxels_tb = np.zeros((length+2, width+2, 1))      
    # adding the generated boundary voxels to the input array
    combined_array_lr = np.concatenate((new_voxels_lr, input_array, new_voxels_lr), axis=0)
    combined_array_lr_bf = np.concatenate((new_voxels_bf, combined_array_lr, new_voxels_bf), axis=1)
    combined_array_lr_bf_tb = np.concatenate((new_voxels_tb, combined_array_lr_bf, new_voxels_tb), axis=2)
    addedboundary_array = combined_array_lr_bf_tb    
    return addedboundary_array

def getkappa_left(calc_w_vec, dx, Tplate_l, Tplate_r):
    ''' calculate the thermal conductivity based on the flow density
    coming from the left side into the resistor network. '''
    length = calc_w_vec.shape[0]
    width = calc_w_vec.shape[1]
    height = calc_w_vec.shape[2]
    
    x_len = (length-2)*dx # length of the structure in m
    Qs = []
    # calculate the average heat flow density from the left side into the resistor network
    for array in calc_w_vec[1, 1:width-1, 1:height-1]:
        for item in array:
            T = item[-1]
            kl = item[2]           
            Q = -kl*(T-Tplate_l)/dx # 1D Fouriers law
            Qs.append(Q)    
    meanQ = np.mean(Qs)          
    # calculate the effective thermal conductivity of the composite
    kappa_sim = -meanQ*(x_len)/(Tplate_r-Tplate_l)
    return kappa_sim

def getresiduals(calc_w_vec):
    ''' calculate the mean residual. '''
    length = calc_w_vec.shape[0]
    width = calc_w_vec.shape[1]
    height = calc_w_vec.shape[2] 
            
    Resid_voxels = []
    # calculate the absolute residual of each node
    for i in range(1, length-1):
        for j in range(1, width-1):
            for k in range(1, height-1):
                
                k_spot,kr,kl,kf,kba,kt,kb,T = calc_w_vec[i, j, k]
                T_right = calc_w_vec[i+1, j, k][-1]
                T_left = calc_w_vec[i-1, j, k][-1]                   
                T_back = calc_w_vec[i, j+1, k][-1]
                T_front = calc_w_vec[i, j-1, k][-1]                 
                T_top = calc_w_vec[i, j, k+1][-1]
                T_bottom = calc_w_vec[i, j, k-1][-1]
                
                Resid = (T*(kr+kl+kf+kba+kt+kb)-T_right*kr-T_left*kl-T_front*kf-T_back*kba-T_top*kt-T_bottom*kb)/(kr+kl+kf+kba+kt+kb) 
                Resid_voxels.append(np.abs(Resid))
    # calcluate the mean residual 
    Mean_Resid = np.mean(Resid_voxels)
    return Mean_Resid  

def kappa_array3D(input_array, k_mat1, k_mat2):    
    ''' input array from composite library, 0's and 1's from the input array 
    correspond to kap_mat1 and kap_mat2 in the resulting array.'''
    kappa_array3D = np.zeros(input_array.shape)
    kappa_array3D[input_array == 0] = k_mat1
    kappa_array3D[input_array == 1] = k_mat2
    return kappa_array3D

def make_calcarray(input_array, Tleft, Tright, Tcomp, dx, R_int):  
    ''' input array: 3D array containing the conductivities of each voxel,
    returns an array containing tuples that store the conductivity of a voxel and
    the thermal conductivities for heat flow in any spacial direction from this spot
    and the actual Temperature of the spot.'''
    length = input_array.shape[0]
    width = input_array.shape[1]
    height = input_array.shape[2]      
    #initialize the array, with tuples containing zeros and T = Tcomp    
    calc_array = np.zeros((length, width, height), dtype=object)    
    for i in range(length):
        for j in range(width): 
            for k in range(height):
                calc_array[i,j,k] = (0,0,0,0,0,0,0,Tcomp)                
    # initialize the hot plate parameters            
    for j in range(1, width-1):
        for k in range(1, height-1):
            # left hot plate
            calc_array[0,j,k] = (input_array[0,j,k],0,0,0,0,0,0,Tleft)
            # right hot plate
            calc_array[-1,j,k] = (input_array[-1,j,k],0,0,0,0,0,0,Tright)           
    # calculate the thermal conductivities for heat flow in each spacial direction for each composite voxel        
    for i in range(1, length-1):
        for j in range(1,width-1):
            for k in range(1,height-1):
                # for each voxel and its neighbours we will now extract the material parameter and actual temperature
                k_spot = input_array[i, j, k]
                k_right = input_array[i+1, j, k]
                k_left = input_array[i-1, j, k]                   
                k_back = input_array[i, j+1, k]
                k_front = input_array[i, j-1, k]                    
                k_top = input_array[i, j, k+1]
                k_bottom = input_array[i, j, k-1] 
                # account for heat flow through different materials and include an interfacial resistance
                if k_right == float('inf'): # exclude interfacial thermal resistances at the boundaries to the hotplates 
                    kr = (0.5*(1/k_spot + 1/k_right))**-1
                elif k_spot == k_right:
                    kr = k_spot                
                elif k_spot != k_right:
                    kr = (0.5*(1/k_spot + 1/k_right) + R_int/dx)**-1
                
                if k_left == float('inf'):
                    kl = (0.5*(1/k_spot + 1/k_left))**-1                
                elif k_spot == k_left:
                    kl = k_spot
                elif k_spot != k_left:
                    kl = (0.5*(1/k_spot + 1/k_left) + R_int/dx)**-1
                
                if k_front == 0: # thermal conductivity at the isolating boundaries is 0
                    kf = 0
                elif k_spot == k_front:
                    kf = k_spot
                elif k_spot != k_front:
                    kf = (0.5*(1/k_spot + 1/k_front) + R_int/dx)**-1
                    
                if k_back == 0:
                    kba = 0    
                elif k_spot == k_back:
                    kba = k_spot
                elif k_spot != k_back:
                    kba = (0.5*(1/k_spot + 1/k_back) + R_int/dx)**-1   
                    
                if k_top == 0:
                    kt = 0
                elif k_spot == k_top:
                    kt = k_spot
                elif k_spot != k_top:
                    kt = (0.5*(1/k_spot + 1/k_top) + R_int/dx)**-1
                    
                if k_bottom == 0:
                    kb = 0
                elif k_spot == k_bottom:
                    kb = k_spot
                elif k_spot != k_bottom:
                    kb = (0.5*(1/k_spot + 1/k_bottom) + R_int/dx)**-1
                    
                calc_array[i,j,k] = (k_spot,kr,kl,kf,kba,kt,kb,Tcomp)            
    return calc_array 

def makelinearplot(x, y, xlabel, ylabel, save=False, name = ''):
    ''' x, y = array,
    xlabel, ylabel, name = str,
    save = bool,
    returns a linear plot of the data.'''
    if len(x) == 1 or len(y) == 1:
        return
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    plt.plot(x,y,'o', color='k', markerfacecolor='salmon', markersize=12, alpha=0.7)      
    plt.xlabel(xlabel, size=15) 
    plt.ylabel(ylabel, size=15)
    plt.yticks(size=12) 
    plt.xticks(size=12) 
    xmin = min(x)-0.05*max(x)
    xmax = max(x)+0.05*max(x) 
    ymin = min(y)-0.05*max(y)
    ymax = max(y)+0.05*max(y) 
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax.set_aspect((xmax-xmin)/(ymax-ymin))  
    ax.tick_params(direction='in', length=7)
    plt.show()
    if save is True:
        fig.savefig(name+'.pdf')

def makelogplot(x, y, xlabel, ylabel, save=False, name=''):
    ''' x, y = array,
    xlabel, ylabel, name = str,
    save = bool,
    returns a linear plot of the data.'''
    if len(x) == 1 or len(y) == 1:
        return
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    plt.plot(x,y,'o', color='k', markerfacecolor='salmon', markersize=12, alpha=0.7)      
    plt.xlabel(xlabel, size=15) 
    plt.ylabel(ylabel, size=15)
    plt.yticks(size=12) 
    plt.xticks(size=12)    
    xmin = min(x)-0.05*max(x)
    xmax = max(x)+0.05*max(x) 
    ymin = 0.1*min(y)
    ymax = 10*max(y)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)  
    ax.set_yscale('log')
    ax.set_aspect((xmax-xmin)/(np.log10(ymax)-np.log10(ymin)))  
    ax.tick_params(direction='in', length=7)
    plt.show()
    if save is True:
        fig.savefig(name+'.pdf')
        
def SOR_update(calc_array, omega=1.979):
    '''updating all Temperatures using the SOR method,
    omega: successive overrelaxation parameter should be choosen between 1 and 2.'''
    length = calc_array.shape[0]
    width = calc_array.shape[1]
    height = calc_array.shape[2]
    
    for i in range(1, length-1):
        for j in range(1, width-1): 
            for k in range(1, height-1):       
                # extract the material parameters and actual temperatures for each voxel and its neighbors
                k_spot,kr,kl,kf,kba,kt,kb,T = calc_array[i, j, k]
                T_right = calc_array[i+1, j, k][-1]
                T_left = calc_array[i-1, j, k][-1]                   
                T_back = calc_array[i, j+1, k][-1]
                T_front = calc_array[i, j-1, k][-1]                 
                T_top = calc_array[i, j, k+1][-1]
                T_bottom = calc_array[i, j, k-1][-1]                              
                # The temperature is updated by solving fouriers law in the spacial directions,
                # so that the amount of heat flow into and out of each node are equal.
                T_new = (T_right*kr+T_left*kl+T_front*kf+T_back*kba+T_top*kt+T_bottom*kb)/(kr+kl+kf+kba+kt+kb)
                # immediately updating the Temperature value in the calc_array                                              
                calc_array[i, j, k] = (k_spot,kr,kl,kf,kba,kt,kb, T + omega*(T_new-T)) 
    return calc_array

def store_x_and_y_data(x, y, x_name, y_name, name=''):
    ''' x,y = array,
    x_name, y_name, name = str, 
    returns txt file containing the data.'''
    header = f"{x_name}\t{y_name}\n"
    data_rows = np.column_stack((x, y))
    allrows = []
    allrows.append(data_rows)
    # Open a dat file to write into
    with open(name+'.dat', 'w') as file:
        # writing the header into the file
        file.write(header)
        # writing the data into the file
        for rows in allrows:
            np.savetxt(file, rows, delimiter='\t')  
    # Close the file
    file.close()

def tempguess_linear(input_array, T_left, T_right):    
    ''' input array is output from the make_calcarray function 
    returns a copy of that array with linearly decreasing initial temperatures from the hot to the cold plate. '''    
    length = input_array.shape[0]
    width = input_array.shape[1]
    height = input_array.shape[2]
    
    calc_linvec_array = np.copy(input_array)
    for i in range(1,length-1):
        for j in range(1,width-1):
            for k in range(1, height-1):
                k_spot,kr,kl,kf,kba,kt,kb,T = input_array[i, j, k]
                # the temperature is dependent on the spacial i-direction
                calc_linvec_array[i, j, k] = (k_spot,kr,kl,kf,kba,kt,kb,(T_left+((T_right-T_left)/(length-2))*(i-1/2))) 
    return calc_linvec_array











             



























 
    
   

    
    

