function S = my_write_output_init(S, filename)

% open .out file
outfname = strcat(filename, '.out'); 
i = 1;
while exist(outfname,'file')
    outfname = sprintf('%s.out_%02d',filename,i);
    i = i + 1;
end

suffixNum = i-1; % save suffix number, only used if suffixNum > 0

% if there are already 100 files, then start using .out only
OUT_MAX = 100;
if i > OUT_MAX
    outfname = strcat(filename,'.out'); 
    suffixNum = -1;
end

% create an output file and write initial variables
fileID = fopen(outfname,'w');
if (fileID == -1) 
    error('\n Cannot open file "%s"\n',outfname);
end 

start_time = fix(clock);
fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'*                      M-SPARC (version Sep 08, 2023)                     *\n');
fprintf(fileID,'*   Copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech   *\n');
fprintf(fileID,'*           Distributed under GNU General Public License 3 (GPL)          *\n');
fprintf(fileID,'*                Date: %s  Start time: %02d:%02d:%02d                  *\n',date,start_time(4),start_time(5),start_time(6));
fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'                           Input parameters                                \n');
fprintf(fileID,'***************************************************************************\n');

if S.Flag_latvec_scale == 0
    fprintf(fileID,'CELL: %f %f %f \n',S.L1,S.L2,S.L3);
    fprintf(fileID,'LATVEC:\n');
    fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_uvec(1,:));
    fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_uvec(2,:));
    fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_uvec(3,:));
else
    fprintf(fileID,'LATVEC_SCALE: %f %f %f \n',S.latvec_scale_x,S.latvec_scale_y,S.latvec_scale_z); 
    fprintf(fileID,'LATVEC:\n');
    fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_vec(1,:));
    fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_vec(2,:));
    fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_vec(3,:));
end

fprintf(fileID,'FD_GRID: %d %d %d\n',S.Nx-S.BCx,S.Ny-S.BCy,S.Nz-S.BCz);
fprintf(fileID,'FD_ORDER: %d\n',S.FDn*2);
%fprintf(fileID,'BOUNDARY_CONDITION: %d\n',S.BC);
str_BC = ['P', 'D'];
fprintf(fileID,'BC:');
fprintf(fileID,' %s',str_BC(S.BCx+1));
fprintf(fileID,' %s',str_BC(S.BCy+1));
fprintf(fileID,' %s',str_BC(S.BCz+1));
fprintf(fileID,'\n');
if (S.BC==2 || S.BC==3 || S.BC==4)
    fprintf(fileID,'KPOINT_GRID: %d %d %d\n',S.nkpt);
    fprintf(fileID,'KPOINT_SHIFT: %d %d %d\n',S.kptshift);
end

if (S.spin_typ ~= 0) 
    fprintf(fileID,'SPIN_TYP: %d\n', S.spin_typ);  
end

if (S.elec_T_type == 0) 
    fprintf(fileID,'ELEC_TEMP_TYPE: fermi-dirac\n');  
elseif (S.elec_T_type == 1) 
    fprintf(fileID,'ELEC_TEMP_TYPE: gaussian\n');  
end
%fprintf(fileID,'ELEC_TEMP: %lf\n',S.elec_T);

fprintf(fileID,'SMEARING: %.9f\n',1/S.bet);
fprintf(fileID,'CHEB_DEGREE: %d\n',S.npl);
fprintf(fileID,'NSTATES: %d\n',S.Nev);
%fprintf(fileID,'NTYPES: %d\n',S.Ntypes);
fprintf(fileID,'EXCHANGE_CORRELATION: %s\n',S.XC);
fprintf(fileID,'CALC_STRESS: %d\n',S.Calc_stress);
if(S.Calc_stress == 0)
    fprintf(fileID,'CALC_PRES: %d\n',S.Calc_pres);
end
%if (S.MDFlag == 1 || S.RelaxFlag == 1)
%    fprintf(fileID,'TWTIME: %f\n',S.TWtime);
%end

if (S.CheFSI_Optmz == 1)
    fprintf(fileID,'CHEFSI_OPTMZ: %d\n',S.CheFSI_Optmz);
end
if (S.chefsibound_flag == 1)
    fprintf(fileID,'CHEFSI_BOUND_FLAG: %d\n',S.chefsibound_flag);
end
if (S.NetCharge ~= 0)
    fprintf(fileID,'NET_CHARGE: %d\n',S.NetCharge);
end
fprintf(fileID,'MAXIT_SCF: %d\n',S.MAXIT_SCF);
if (S.MDFlag == 1)
    fprintf(fileID,'MD_FLAG: %d\n',S.MDFlag);
    fprintf(fileID,'MD_METHOD: %s\n',S.MDMeth);
    fprintf(fileID,'MD_TIMESTEP: %.2f\n',S.MD_dt); 
    %fprintf(fileID,'ATOMIC_MASS:');
    %for (i = 0; i < S.Ntypes; i++)     
    %     fprintf(fileID,' %.15f', S.Mass[i]);      
    %end
    fprintf(fileID,'MD_NSTEP: %d\n',S.MD_Nstep);
    fprintf(fileID,'ION_ELEC_EQT: %d\n',S.ion_elec_eqT);
    % fprintf(fileID,'ION_VEL_DSTR: %d\n',S.ion_vel_dstr);
    fprintf(fileID,'ION_TEMP: %f\n',S.ion_T);
    % if(strcmp(S.MDMeth,'NVT_NH'))
        % fprintf(fileID,'ION_TEMP_END: %lf\n',S.thermos_Tf);
        % fprintf(fileID,'QMASS: %lf\n',S.qmass);
    % end
end
if (S.RelaxFlag==1)
    fprintf(fileID,'RELAX_FLAG: %d\n',S.RelaxFlag);
    fprintf(fileID,'RELAX_METHOD: %s\n',S.RelaxMeth);
    fprintf(fileID,'RELAX_NITER: %d\n',S.max_relax_it);
    if(strcmp(S.RelaxMeth,'LBFGS'))
        fprintf(fileID,'L_HISTORY: %d\n',S.L_history);
        fprintf(fileID,'L_FINIT_STP: %f\n',S.L_finit_stp);
        fprintf(fileID,'L_MAXMOV: %f\n',S.L_maxmov);
        fprintf(fileID,'L_AUTOSCALE: %d\n',S.L_autoscale);
        fprintf(fileID,'L_LINEOPT: %d\n',S.L_lineopt);
        fprintf(fileID,'L_ICURV: %f\n',S.L_icurv);
    elseif (strcmp(S.RelaxMeth,'NLCG'))
        fprintf(fileID,'NLCG_SIGMA: %f\n',S.NLCG_sigma);
    elseif (strcmp(S.RelaxMeth,'FIRE'))
        fprintf(fileID,'FIRE_dt: %f\n',S.FIRE_dt);
        fprintf(fileID,'FIRE_mass: %f\n',S.FIRE_mass);
        fprintf(fileID,'FIRE_maxmov: %f\n',S.FIRE_maxmov);
    end
    fprintf(fileID,'TOL_RELAX: %.2E\n',S.TOL_RELAX);
elseif (S.RelaxFlag==2)    
    fprintf(fileID,'RELAX_FLAG: %d\n',S.RelaxFlag);    
    fprintf(fileID,'RELAX_NITER: %d\n',S.max_relax_it);    
    fprintf(fileID,'TOL_RELAX_CELL: %.2E\n',S.TOL_RELAX_CELL);    
    fprintf(fileID,'RELAX_MAXDILAT: %f\n',S.max_dilatation);    
end
fprintf(fileID,'TOL_SCF: %.2E\n',S.SCF_tol);
fprintf(fileID,'TOL_POISSON: %.2E\n',S.poisson_tol);
fprintf(fileID,'TOL_LANCZOS: %.2E\n',S.TOL_LANCZOS);
fprintf(fileID,'TOL_PSEUDOCHARGE: %.2E\n',S.pseudocharge_tol);
if (S.MixingVariable == 0)
    fprintf(fileID,'MIXING_VARIABLE: density\n');
elseif  (S.MixingVariable == 1)
    fprintf(fileID,'MIXING_VARIABLE: potential\n');
end
if (S.MixingPrecond == 0)
    fprintf(fileID,'MIXING_PRECOND: none\n');
elseif (S.MixingPrecond == 1)
    fprintf(fileID,'MIXING_PRECOND: kerker\n');
% elseif (S.MixingPrecond == 2)
%     fprintf(fileID,'MIXING_PRECOND: resta\n');
% elseif (S.MixingPrecond == 3)
%     fprintf(fileID,'MIXING_PRECOND: truncated_kerker\n');
end

% for large periodic systems, give warning if preconditioner is not chosen
if S.BC == 2 || 0
    L_diag = sqrt(S.L1^2 + S.L2^2 + S.L3^2);
    if L_diag > 20 && S.MixingPrecond == 0
        fprintf(fileID,"#WARNING: the preconditioner for SCF has been turned off, this \n");
        fprintf(fileID,"might lead to slow SCF convergence. To specify SCF preconditioner, \n");
        fprintf(fileID,"#use 'MIXING_PRECOND' in the .inpt file\n");
    end
end
if S.spin_typ ~= 0
    if (S.MixingPrecondMag == 0)
        fprintf(fileID,'MIXING_PRECOND_MAG: none\n');
    elseif (S.MixingPrecondMag == 1)
        fprintf(fileID,'MIXING_PRECOND_MAG: kerker\n');
%     elseif (S.MixingPrecondMag == 2)
%         fprintf(fileID,'MIXING_PRECOND_MAG: resta\n');
%     elseif (S.MixingPrecondMag == 3)
%         fprintf(fileID,'MIXING_PRECOND_MAG: truncated_kerker\n');
    end
end
if (S.MixingPrecond ~= 0)
    fprintf(fileID,'TOL_PRECOND: %.2E\n',S.precond_tol);
end
if (S.MixingPrecond == 1) % kerker
    fprintf(fileID,'PRECOND_KERKER_KTF: %.2f\n',S.precond_kerker_kTF);
    fprintf(fileID,'PRECOND_KERKER_THRESH: %.2f\n',S.precond_kerker_thresh);
% elseif (S.MixingPrecond == 2) % resta
%     %fprintf(fileID,'TOL_PRECOND: %.2E\n',S.precond_tol);
%     fprintf(fileID,'PRECOND_RESTA_Q0: %.3f\n',S.precond_resta_q0);
%     fprintf(fileID,'PRECOND_RESTA_RS: %.3f\n',S.precond_resta_Rs);
%     fprintf(fileID,'PRECOND_FITPOW: %d\n',S.precond_fitpow);
% elseif (S.MixingPrecond == 3) % truncated kerker
%     fprintf(fileID,'PRECOND_KERKER_KTF: %.2f\n',S.precond_kerker_kTF);
%     fprintf(fileID,'PRECOND_KERKER_THRESH: %.2f\n',S.precond_kerker_thresh);
%     fprintf(fileID,'PRECOND_FITPOW: %d\n',S.precond_fitpow);
end
if S.spin_typ ~= 0
    if S.MixingPrecondMag == 1
        fprintf(fileID,'PRECOND_KERKER_KTF_MAG: %.2f\n',S.precond_kerker_kTF_mag);
        fprintf(fileID,'PRECOND_KERKER_THRESH_MAG: %.2f\n',S.precond_kerker_thresh_mag);
    end
end
fprintf(fileID,'MIXING_PARAMETER: %.2f\n',S.MixingParameter);
if S.PulayFrequency > 1
    fprintf(fileID,'MIXING_PARAMETER_SIMPLE: %.2f\n',S.MixingParameterSimple);
end
if S.spin_typ ~= 0
    fprintf(fileID,'MIXING_PARAMETER_MAG: %.2f\n',S.MixingParameterMag);
    if S.PulayFrequency > 1
        fprintf(fileID,'MIXING_PARAMETER_SIMPLE_MAG: %.2f\n',S.MixingParameterSimpleMag);
    end
end
fprintf(fileID,'MIXING_HISTORY: %d\n',S.MixingHistory);
fprintf(fileID,'PULAY_FREQUENCY: %d\n',S.PulayFrequency);
fprintf(fileID,'PULAY_RESTART: %d\n',S.PulayRestartFlag);
fprintf(fileID,'REFERENCE_CUTOFF: %.2f\n',S.rc_ref);
fprintf(fileID,'RHO_TRIGGER: %d\n',S.rhoTrigger);
fprintf(fileID,'NUM_CHEFSI: %d\n',S.nchefsi);
fprintf(fileID,'FIX_RAND: %d\n',S.FixRandSeed);
fprintf(fileID,'PRINT_FORCES: %d\n',S.PrintForceFlag);
fprintf(fileID,'PRINT_ATOMS: %d\n',S.PrintAtomPosFlag);
fprintf(fileID,'PRINT_EIGEN: %d\n',S.PrintEigenFlag);
fprintf(fileID,'PRINT_DENSITY: %d\n',S.PrintElecDensFlag);
if(S.MDFlag == 1)
    fprintf(fileID,'PRINT_MDOUT: %d\n',S.PrintMDout);
end
if(S.MDFlag == 1 || S.RelaxFlag == 1)  
    fprintf(fileID,'PRINT_VELS: %d\n',S.PrintAtomVelFlag);  
    fprintf(fileID,'PRINT_RESTART: %d\n',S.Printrestart);
    if(S.Printrestart == 1)
        fprintf(fileID,'PRINT_RESTART_FQ: %d\n',S.Printrestart_fq);
    end
end

if(S.RelaxFlag == 1)
    fprintf(fileID,'PRINT_RELAXOUT: %d\n',S.PrintRelaxout);
end

fprintf(fileID,'OUTPUT_FILE: %s\n',outfname);
if (S.RestartFlag == 1)
    fprintf(fileID,'RESTART_FLAG: %d\n',S.RestartFlag);
end

if(S.usefock == 1)
    fprintf(fileID,'MAXIT_FOCK: %d\n',S.MAXIT_FOCK);
    fprintf(fileID,'MINIT_FOCK: %d\n',S.MINIT_FOCK);
    fprintf(fileID,'TOL_FOCK: %.2E\n',S.FOCK_TOL);
    fprintf(fileID,'TOL_SCF_INIT: %.2E\n',S.SCF_tol_init);
    fprintf(fileID,'EXX_METHOD: %s\n',S.ExxMethod);
    fprintf(fileID,'ACE_FLAG: %d\n',S.ACEFlag);
    if S.ACEFlag == 1
        fprintf(fileID,'EXX_ACE_VALENCE_STATES: %d\n',S.EXXACEVal_state);
    end
    if S.BC == 2
        fprintf(fileID,'EXX_DOWNSAMPLING: %d %d %d\n',S.exx_downsampling);
    end
    fprintf(fileID,'EXX_DIVERGENCE: %s\n', S.ExxDivMethod);
    if S.xc == 427
        fprintf(fileID,'EXX_RANGE_FOCK: %.4f\n', S.hyb_range_fock);
        fprintf(fileID,'EXX_RANGE_PBE: %.4f\n', S.hyb_range_pbe);
    end
end

if(S.d3Flag == 1)
    fprintf(fileID,'D3_FLAG: %d\n',S.d3Flag);
    fprintf(fileID,'D3_RTHR: %f\n',S.d3Rthr);
    fprintf(fileID,'D3_CN_THR: %f\n',S.d3Cn_thr);
end

if(S.vdWDFFlag == 1 || S.vdWDFFlag == 2)  
    fprintf(fileID,'VDWDF_GEN_KERNEL: %d\n',S.vdWDFKernelGenFlag);  
end

fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'                                Cell                                       \n');
fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'Lattice vectors:\n');
fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_uvec(1,:)*S.L1);
fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_uvec(2,:)*S.L2);
fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_uvec(3,:)*S.L3);
fprintf(fileID,'Volume :%18.10E (Bohr^3)\n',S.L1*S.L2*S.L3*S.Jacb);
fprintf(fileID,'Density :%18.10E (amu/Bohr^3), %18.10E (g/cc)\n',...
            S.TotalMass/(S.L1*S.L2*S.L3*S.Jacb), S.TotalMass/(S.L1*S.L2*S.L3*S.Jacb)*11.2058730627683);

% fprintf(fileID,'***************************************************************************\n');
% fprintf(fileID,'                           Parallelization                                 \n');
% fprintf(fileID,'***************************************************************************\n');
% fprintf(fileID,'NP_KPOINT_PARAL: %d\n',S.npkpt);
% fprintf(fileID,'NP_BAND_PARAL: %d\n',S.npband);
% fprintf(fileID,'NP_DOMAIN_PARAL: %d %d %d\n',S.npNdx,S.npNdy,S.npNdz);
% fprintf(fileID,'NP_DOMAIN_PHI_PARAL: %d %d %d\n',S.npNdx_phi,S.npNdy_phi,S.npNdz_phi);

fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'                             Initialization                                \n');
fprintf(fileID,'***************************************************************************\n');
% fprintf(fileID,'Number of processors               :  %d\n',nproc);

if ( (abs(S.dx-S.dy) <=1e-12) && (abs(S.dx-S.dz) <=1e-12) ...
    && (abs(S.dy-S.dz) <=1e-12) ) 
    fprintf(fileID,'Mesh spacing                       : % f (Bohr)\n',S.dx); 
else
    fprintf(fileID,'Mesh spacing in x-direction        : % f (Bohr)\n',S.dx); 
    fprintf(fileID,'Mesh spacing in y-direction        : % f (Bohr)\n',S.dy); 
    fprintf(fileID,'Mesh spacing in z direction        : % f (Bohr)\n',S.dz);     
end

if (S.BC==2 || S.BC==3 || S.BC==4) 
    fprintf(fileID,'Number of symmetry adapted k-points:  %d\n',S.tnkpt);       
end

fprintf(fileID,'Output printed to                  :  %s\n',outfname);

%if (S.PrintAtomPosFlag==1)
%    fprintf(fileID,'Atom positions printed to          :  %s\n',S.AtomFilename);      
%end

%if (S.PrintForceFlag==1)
%    fprintf(fileID,'Forces printed to                  :  %s\n',S.ForceFilename);
%end 

% if (S.PrintEigenFlag==1)
%     fprintf(fileID,'Final eigenvalues printed to       :  %s\n',S.EigenFilename);
% end

% if (S.MDFlag == 1 && S.PrintMDout == 1)
%     fprintf(fileID,'MD output printed to               :  %s\n',S.MDFilename);
% end

% if (S.RelaxFlag == 1 && S.PrintRelaxout == 1)
%     fprintf(fileID,'Relax output printed to            :  %s\n',S.RelaxFilename);
% end

fprintf(fileID,'Total number of atom types         :  %d\n',S.n_typ);
fprintf(fileID,'Total number of atoms              :  %d\n',S.n_atm);
fprintf(fileID,'Total number of electrons          :  %d\n',S.Nelectron);

for ityp = 1:S.n_typ
    fprintf(fileID,'Atom type %-2d (valence electrons)   :  %s %d\n',ityp,S.Atm(ityp).typ, S.Atm(ityp).Z);
    fprintf(fileID,'Pseudopotential                    :  %s\n',S.Atm(ityp).psdfname);     
    fprintf(fileID,'lloc                               :  %d\n',S.Atm(ityp).lloc);    
    fprintf(fileID,'Atomic mass                        :  %.15f\n',S.Atm(ityp).Mass);
    fprintf(fileID,'Pseudocharge radii of atom type %-2d :  %.2f %.2f %.2f\n',ityp,S.Atm(ityp).rb_x,S.Atm(ityp).rb_y,S.Atm(ityp).rb_z);
    fprintf(fileID,'Number of atoms of type %-2d         :  %d\n',ityp,S.Atm(ityp).n_atm_typ);
    % if (S.PrintAtomPosFlag == 1 && S.MDFlag == 0 && S.RelaxFlag == 0)
    %     fprintf(fileID,'Fractional coordinates of atoms of type %-2d    :\n',ityp);
    %     for j = 1:S.Atm(ityp).n_atm_typ
    %         fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atm(ityp).coords(j,:)./[S.L1, S.L2, S.L3]);
    %     end
    % end
end

[~, mem_num, mem_unit] = my_print_mem(S.memory_usage);
fprintf(fileID, 'Estimated total memory usage       :  %-.2f %s\n', mem_num, mem_unit);

fclose(fileID);    

%-----------------------------------------
% Write atom positions to .static file 
%-----------------------------------------
if ((S.PrintAtomPosFlag == 1 || S.PrintForceFlag == 1) && S.MDFlag == 0 && S.RelaxFlag == 0)
    staticfname = strcat(filename,'.static'); 
    if suffixNum > 0
        staticfname = sprintf('%s.static_%02d',filename,suffixNum);
    end
    % open file
    fid = fopen(staticfname,'w') ;
    assert(fid~=-1,'Error: Cannot open .static file %s',staticfname);

    if(S.PrintAtomPosFlag == 1)
        fprintf(fid,'***************************************************************************\n');
        fprintf(fid,'                            Atom positions                                 \n');
        fprintf(fid,'***************************************************************************\n');

        nFracCoord = sum(S.IsFrac);
        
        if nFracCoord == S.n_typ
            fprintf(fileID,'Fractional coordinates of atoms:\n');
            for j = 1:S.n_atm
                fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atoms(j,:)./[S.L1, S.L2, S.L3]);
            end
        elseif nFracCoord == 0
            fprintf(fileID,'Cartesian coordinates of atoms (Bohr):\n');
            for j = 1:S.n_atm
                atomJPos = S.Atoms(j,:)*S.lat_uvec;
                fprintf(fileID,'%18.10f %18.10f %18.10f\n',atomJPos(1), atomJPos(2), atomJPos(3));
            end
        else
            for ityp = 1:S.n_typ
                if S.IsFrac(ityp)
                    fprintf(fileID,'Fractional coordinates of %s:\n',S.Atm(ityp).typ); 
                    for j = 1:S.Atm(ityp).n_atm_typ
                        fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atm(ityp).coords(j,:)./[S.L1, S.L2, S.L3]);
                    end
                else
                    fprintf(fileID,'Cartesian coordinates of %s (Bohr):\n',S.Atm(ityp).typ); 
                    for j = 1:S.Atm(ityp).n_atm_typ
                        atomJPos = S.Atm(ityp).coords(j,:)*S.lat_uvec;
                        fprintf(fileID,'%18.10f %18.10f %18.10f\n',atomJPos(1), atomJPos(2), atomJPos(3));
                    end
                end
            end
        end
    end

    if S.spin_typ ~= 0
        fprintf(fileID,'Initial spin:\n');
        for ityp = 1:S.n_typ
            for j = 1:S.Atm(ityp).n_atm_typ
                fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atm(ityp).mag(j,:));
            end
        end
    end

    % close file
    fclose(fid);
    
    S.staticfname = staticfname; % save to structure;
end


% save to structure
S.suffixNum   = suffixNum;
S.outfname    = outfname; 

end
