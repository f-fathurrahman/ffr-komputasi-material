clear all; close all

t = 0.5;
U = 0.1;
noOfSites = 2;
noOfUp = 1;
noOfDn = 1;

% Generate basis
[ combinedBasis, totalNoOfPossiblestates,totalNoOfUpStates, ...
  totalNoOfDnStates, upStates, dnStates ] = ...
  generateBasis( noOfSites, noOfUp, noOfDn );

max_kinetic_num_non_zero_elems = 2*(noOfSites - noOfUp - 2) + 2*(noOfSites - noOfDn - 2);
fprintf('max_kinetic_num_non_zero_elems = %d\n', max_kinetic_num_non_zero_elems)
KINETIC_COUNTER = 0;
kinetic_rows = zeros(max_kinetic_num_non_zero_elems, 1);
kinetic_cols = zeros(max_kinetic_num_non_zero_elems, 1);
kinetic_elems = zeros(max_kinetic_num_non_zero_elems, 1);

% the number of electrons to be hopped over if the electrons hop
% around the lattice boundary (can easily see that this must be the case):
noOfUpInterior = noOfUp - 1;
noOfDnInterior = noOfDn - 1;

for m=1:totalNoOfPossiblestates % go through each state in the basis:

  % save the unshifted spin up and spin down sectors:
  upSectorDec = combinedBasis(m, 2);
  dnSectorDec = combinedBasis(m, 3);
  
  upSector= de2bi_modified(upSectorDec, noOfSites);
  dnSector= de2bi_modified(dnSectorDec, noOfSites);
  
  % find the occupied lattice sites:    
  upNonZero = find(upSector);
  dnNonZero = find(dnSector);

  fprintf('\n--------------------------\n')
  fprintf('m = %d\n', m)
  disp('upSector'); disp(upSector)
  disp('dnSector'); disp(dnSector)
  disp('combinedBasis'); disp(combinedBasis(m,2:3))
  fprintf('--------------------------\n')

  % shift for spin up:
  for n = upNonZero % for each occupied site

    % left shift:
    leftShiftResult = upSector;
    % figure out which site is the one to its left (periodic boundary condition)
    leftShiftedIndex = mod(n-2,noOfSites) + 1;

    fprintf('\nSPIN UP m = %3d n = %3d\n', m, n);

    if upSector(leftShiftedIndex) ~= 1

      % perform the shift:
      leftShiftResult(n) = 0;
      leftShiftResult(leftShiftedIndex) = 1;

      disp('Left shift is done');
      disp('leftShiftResult = ');
      disp(leftShiftResult);
           
      % figure out where in the basis this shifted state is
      upIndexOfLeftShiftedResult = binaraysearchasc(upStates, bi2de_modified(leftShiftResult) );
      dnIndexOfLeftShiftedResult = mod(m-1,totalNoOfDnStates) + 1;
      basisIndexOfLeftShiftedResult = (upIndexOfLeftShiftedResult-1)*totalNoOfDnStates + ...
          dnIndexOfLeftShiftedResult;
        
      fprintf('leftShiftedIndex=%d n=%d\n', leftShiftedIndex, n)
      % update that state:
      if leftShiftedIndex < n % if the electron does not hop around the boundary
        KINETIC_COUNTER = KINETIC_COUNTER + 1;
        kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
        kinetic_cols(KINETIC_COUNTER) = m;
        kinetic_elems(KINETIC_COUNTER) = -t;
        fprintf('Does not hop the doundary: KINETIC_COUNTER = %d\n', KINETIC_COUNTER)             
      else % if the electron does hop around the boundary
        if mod(noOfUpInterior,2)== 0 % if the number of electrons to be hopped over is even
          KINETIC_COUNTER = KINETIC_COUNTER + 1;
          kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
          kinetic_cols(KINETIC_COUNTER) = m;
          kinetic_elems(KINETIC_COUNTER) = -t;
        else
          KINETIC_COUNTER = KINETIC_COUNTER + 1;
          kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
          kinetic_cols(KINETIC_COUNTER) = m;
          kinetic_elems(KINETIC_COUNTER) = +t;
        end
        fprintf('Hops the doundary: KINETIC_COUNTER = %d\n', KINETIC_COUNTER)
      end 
    else
      disp('Left shift is not done')
    end
       
    % right shift:
    rightShiftResult = upSector;
    rightShiftedIndex = mod(n,noOfSites) + 1;
    if upSector(rightShiftedIndex) ~= 1
      % perform the shift:
      rightShiftResult(n)=0;           
      rightShiftResult(rightShiftedIndex)=1;

      disp('Right shift is done');
      disp('rightShiftResult = ');
      disp(rightShiftResult);

      % figure out where in the basis this shifted state is
      upIndexOfRightShiftedResult = binaraysearchasc(upStates, bi2de_modified(rightShiftResult) );
      dnIndexOfRightShiftedResult = mod(m-1,totalNoOfDnStates) + 1;
      basisIndexOfRightShiftedResult = (upIndexOfRightShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;
  
      fprintf('rightShiftedIndex=%d n=%d\n', rightShiftedIndex, n)
      % update that state:
      if rightShiftedIndex > n
        KINETIC_COUNTER = KINETIC_COUNTER + 1;
        kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
        kinetic_cols(KINETIC_COUNTER) = m;
        kinetic_elems(KINETIC_COUNTER) = -t;
        fprintf('Hops the doundary: KINETIC_COUNTER = %d\n', KINETIC_COUNTER)
      else
        if mod(noOfUpInterior,2)== 0 
          KINETIC_COUNTER = KINETIC_COUNTER + 1;
          kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
          kinetic_cols(KINETIC_COUNTER) = m;
          kinetic_elems(KINETIC_COUNTER) = -t;
        else
          KINETIC_COUNTER = KINETIC_COUNTER + 1;
          kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
          kinetic_cols(KINETIC_COUNTER) = m;
          kinetic_elems(KINETIC_COUNTER) = +t;
        end
        fprintf('Does not hop the doundary: KINETIC_COUNTER = %d\n', KINETIC_COUNTER)
      end
    else
      disp('Right shift is not done')    
    end
  
  end % for  

  % shift for spin down:
  for p = dnNonZero

    fprintf('\nSPIN DN m = %3d p = %3d\n', m, p);

    % left shift:
    leftShiftResult=dnSector;
    leftShiftedIndex=mod( p-2,noOfSites) + 1; 
    if dnSector(leftShiftedIndex) ~= 1
      % perform the shift:
      leftShiftResult(p)=0;           
      leftShiftResult(leftShiftedIndex)=1;

      disp('Left shift is done');
      disp('leftShiftResult = ');
      disp(leftShiftResult);

      % figure out where in the basis this shifted state is
      dnIndexOfLeftShiftedResult =  binaraysearchasc(dnStates, bi2de_modified(leftShiftResult) );
      upIndexOfLeftShiftedResult = floor(( m - 1 )/totalNoOfDnStates)+1;
      basisIndexOfLeftShiftedResult = (upIndexOfLeftShiftedResult-1)*totalNoOfDnStates + ...
        dnIndexOfLeftShiftedResult;
           
      if leftShiftedIndex < p % if the electron does not hop around the boundary
        disp('Electron does not hop the doundary')
        KINETIC_COUNTER = KINETIC_COUNTER + 1;
        kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
        kinetic_cols(KINETIC_COUNTER) = m;
        kinetic_elems(KINETIC_COUNTER) = -t;
      else % if the electron does hop around the boundary      
        disp('Electron hops the doundary')         
        if mod(noOfDnInterior,2)== 0 % if that number is even
            KINETIC_COUNTER = KINETIC_COUNTER + 1;
            kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
            kinetic_cols(KINETIC_COUNTER) = m;
            kinetic_elems(KINETIC_COUNTER) = -t;
        else
            KINETIC_COUNTER = KINETIC_COUNTER + 1;
            kinetic_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
            kinetic_cols(KINETIC_COUNTER) = m;
            kinetic_elems(KINETIC_COUNTER) = +t;
        end
      end
    end
      
    % right shift:
    rightShiftResult = dnSector;
    rightShiftedIndex = mod(p,noOfSites)+1;  
       
    if dnSector(rightShiftedIndex) ~= 1
      % perform the shift:
      rightShiftResult(p) = 0;
      rightShiftResult(rightShiftedIndex)=1;
      
      disp('Right shift is done');
      disp('rightShiftResult = ');
      disp(rightShiftResult);

      % figure out where in the basis this shifted state is
      dnIndexOfRightShiftedResult = binaraysearchasc(dnStates, bi2de_modified(rightShiftResult) );
      upIndexOfLeftShiftedResult = floor(( m - 1 )/totalNoOfDnStates)+1;
      basisIndexOfRightShiftedResult = (upIndexOfLeftShiftedResult-1)*totalNoOfDnStates + ...
        dnIndexOfRightShiftedResult;
         
      if rightShiftedIndex > p % if the electron does not hop around the boundary
        disp('Electron does not hop the doundary')
        KINETIC_COUNTER = KINETIC_COUNTER + 1;
        kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
        kinetic_cols(KINETIC_COUNTER) = m;
        kinetic_elems(KINETIC_COUNTER) = -t;
      else % if the electron does hop around the boundary
        disp('Electron hops the doundary')
        if mod(noOfDnInterior,2)== 0 % if that number is even
          KINETIC_COUNTER = KINETIC_COUNTER + 1;
          kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
          kinetic_cols(KINETIC_COUNTER) = m;
          kinetic_elems(KINETIC_COUNTER) = -t;
        else
          KINETIC_COUNTER = KINETIC_COUNTER + 1;
          kinetic_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
          kinetic_cols(KINETIC_COUNTER) = m;
          kinetic_elems(KINETIC_COUNTER) = t;
        end
      end
    end
  end    
end

kinetic_rows = kinetic_rows( kinetic_rows ~= 0);
kinetic_cols = kinetic_cols( 1:length(kinetic_rows));
kinetic_elems = kinetic_elems( 1:length(kinetic_rows));
kineticHamiltonian = sparse( kinetic_rows, kinetic_cols, ...
  kinetic_elems, totalNoOfPossiblestates, totalNoOfPossiblestates);
