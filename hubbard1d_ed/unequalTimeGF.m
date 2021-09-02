function [ spinUpGreenFunction, spinDnGreenFunction ] = ...
  unequalTimeGF( t, U, tau, noOfSites, noOfUp, noOfDn )

  % calculate the equal time GF by constructing separate matrices
  % for the c_i and c_j^\dagger operators

  tic;

  do_it = (noOfUp < noOfSites) && (noOfDn < noOfSites);
  if ~do_it
    error('Error: cannot apply creation operator when number of electrons = number of sites');
  end

  spinUpGreenFunction=zeros(noOfSites);
  spinDnGreenFunction=zeros(noOfSites);

  %ASSUMING THAT THE HAMILTONIAN IS REAL SYMMETRIC
  [groundState, groundStateEnergy] = ...
    eigs( hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn ), 1,'sa'); 
       
  %% SPIN UP:
  sizeSpacePlusOne = nchoosek(noOfSites,noOfUp+1)*nchoosek(noOfSites,noOfDn);
  % the Hamiltonian in expanded space:
  expmSecondHamiltonian = expm( -tau * speye(sizeSpacePlusOne) * ...
                                 hubbardHamiltonian( t, U, noOfSites, noOfUp+1, noOfDn ) );
  for i=1:noOfSites        
    % destruction operator:
    % need to take Hermitian conjugate to turn creation into destruction operator
    destructionMatrix = creationOperator( noOfSites, noOfUp, noOfDn , i, 'up' )'; 
            
    for j=1:noOfSites
      right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j, 'up' ) * ...
                                              groundState;
      left_wave_function =  (groundState') * destructionMatrix; 
      spinUpGreenFunction(i,j) = (left_wave_function * expmSecondHamiltonian * right_wave_function) * exp(tau*groundStateEnergy);
            
      clearvars left_wave_function right_wave_function;
    end
        
    clearvars destructionMatrix;
  end
    
  clearvars secondHamiltonian expmSecondHamiltonian;
    
  %% SPIN DOWN:
  sizeSpacePlusOne = nchoosek(noOfSites,noOfUp)*nchoosek(noOfSites,noOfDn+1);  
  % the Hamiltonian in expanded space:
  expmSecondHamiltonian = expm( -tau * speye(sizeSpacePlusOne) * ...
                                   hubbardHamiltonian( t, U, noOfSites, noOfUp, noOfDn+1 ));
  for i=1:noOfSites        
    % destruction operator:
    % need to take Hermitian conjugate to turn creation into destruction operator
    destructionMatrix = creationOperator( noOfSites, noOfUp, noOfDn , i, 'dn' )'; 
            
    for j=1:noOfSites
      right_wave_function = creationOperator( noOfSites, noOfUp, noOfDn , j, 'dn' ) * ...
                                              groundState;
            
      left_wave_function =  (groundState') * destructionMatrix;
      spinDnGreenFunction(i,j) = (left_wave_function * expmSecondHamiltonian * right_wave_function) * exp(tau*groundStateEnergy); 
            
      clearvars left_wave_function right_wave_function;
    end
        
    clearvars destructionMatrix;
  end    
    
  clearvars secondHamiltonian expmSecondHamiltonian;

  elapsed_time = toc;

  fprintf('Elapsed time = %f\n', elapsed_time)

end
