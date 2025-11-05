import os
import composite_lib as comp
import hotplate_lib as hot


'''input parameters to run the calculation'''
name = input('samplename: ')
Lvox = 60 # length of total structure in voxels
disp_volfrac = int(input('volumefraction of the dispersed phase / %: '))
sclust = 4 # number of voxels per cluster of the dispersed Phase

kap_cont = 0.71 # W/(m*K), thermal cond (continuous Phase)  
kap_disp = 0.32 # W/(m*K), thermal cond (dispersed Phase)  
R_int = 1/(500000) # K*m^2/W, interfacial thermal resistance

x_len = 600*10**-6 # m, total length of the structure
dx = x_len/Lvox # m, distance between two voxels in m

T_left = 50 # °C, left hotplate Temperature
T_right = 0 # °C, right hotplate Temperature
T_comp = 25 # °C, initial composite Temperature

###############################################################################

'''generating and plotting the inital voxel structure''' 
voxelstructure = comp.blobs3D(Lvox, disp_volfrac, sclust) # generate the voxelstructure 
comp.compositefigure(voxelstructure, show=True, save=True, name=name) # plot and save the figure of the voxelstructure

###############################################################################

'''generating the RN'''
kappa_array = hot.kappa_array3D(voxelstructure, kap_cont, kap_disp) # assigning conductivities to the voxels representing different phases
kappa_array_w_bounds = hot.add_boundaries(kappa_array) # building hot plates (left and right) and adiabatic conditions (top, bottom, front, back)
calc_array = hot.make_calcarray(kappa_array_w_bounds, T_left, T_right, T_comp, dx, R_int) # building the resistor network
calc_array_w_vec = hot.tempguess_linear(calc_array, T_left, T_right) # initializing the Temperatures with temperatures linearly decreasing from the hot to the cold plate

###############################################################################

'''approximating the steady state'''
max_iter = 100000 
cutoff_resid = 1*10**-9 

iteration = 0
iterations = []
residuals = []
kappas_sim = []
   
# iterating using the SOR method until cutoff criterion or max number of iterations is reached
while iteration < max_iter:    
    # calculating residual and kappa of the initital temperature distribution    
    if iteration == 0:
        iterations.append(iteration)        
        Mean_Resid = hot.getresiduals(calc_array_w_vec)         
        residuals.append(Mean_Resid)        
        kappa_sim_left = hot.getkappa_left(calc_array_w_vec, dx, T_left, T_right)    
        kappas_sim.append(kappa_sim_left)
        # Save data, in case the initial condition already fulfils the cuttoff criteria
        if residuals[-1] < cutoff_resid:
            hot.store_x_and_y_data(iterations, kappas_sim, 'iteration','kappa_sim', name=f'{name}_kappas_sim_finalkappa_{kappa_sim_left:.5f}')
            hot.store_x_and_y_data(iterations, residuals, 'iteration', 'residual', name=f'{name}_residuals_finalresid_{Mean_Resid:.2e}')
            print('iteration:', iteration) 
            print('k_sim:', f'{kappa_sim_left:.4f}')
            print('Mean Residual:', f'{Mean_Resid:.2e}') 
            break      
       
    iteration += 1     
    calc_array_w_vec = hot.SOR_update(calc_array_w_vec, omega=1.979) # updating all temperatures using the SOR method  
    
    if iteration % 10 == 0:
        print('iteration:', iteration)   
    if iteration % 50 == 0:
        iterations.append(iteration)
        # calculate kappa and make a plot
        kappa_sim_left = hot.getkappa_left(calc_array_w_vec, dx, T_left, T_right)        
        kappas_sim.append(kappa_sim_left)
        print('k_sim:', f'{kappa_sim_left:.4f}')
        hot.makelinearplot(iterations, kappas_sim,'iteration', r'$\kappa \,/\, W \cdot (\text{m} \cdot \text{K})^{-1}$', save=False)
        # calculate residual and make a plot
        Mean_Resid = hot.getresiduals(calc_array_w_vec)      
        residuals.append(Mean_Resid)
        print('Mean Residual:', f'{Mean_Resid:.2e}') 
        hot.makelogplot(iterations, residuals, 'iteration', 'residual / K', save=False)
        
        # this can help at the HPC
        curr_direc = os.getcwd()
        for filename in os.listdir(curr_direc):
            if '_act' in filename and name in filename:
                file_path = os.path.join(curr_direc, filename)
                os.remove(file_path)
        if residuals[-1] > cutoff_resid:
            hot.store_x_and_y_data(iterations, kappas_sim, 'iteration','kappa_sim', name=f'{name}_iter_{iteration}_kappas_sim_actkappa_{kappa_sim_left:.5f}')
            hot.store_x_and_y_data(iterations, residuals, 'iteration', 'residual', name=f'{name}_residuals_actresid_{Mean_Resid:.2e}')     
        
    # save data if maximum number of iterations is reached or the cuttoff criterion is met
    if residuals[-1] < cutoff_resid or iteration == max_iter:
        hot.store_x_and_y_data(iterations, kappas_sim, 'iteration','kappa_sim', name=f'{name}_kappas_sim_finalkappa_{kappa_sim_left:.5f}')
        hot.store_x_and_y_data(iterations, residuals, 'iteration', 'residual', name=f'{name}_residuals_finalresid_{Mean_Resid:.2e}')
        hot.makelogplot(iterations, residuals, 'iteration', 'residual / K', save=True, name=name+'_residualplot')
        break  



