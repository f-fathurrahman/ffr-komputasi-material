function S = my_check_atomlocation(S)

% map atom positions back to the domain if atoms are outside the domain in
% periodic directions or throw an error in Dirichlet directions
if S.BCx == 1
    if (sum(S.Atoms(:,1) >= S.L1 | S.Atoms(:,1) < 0) > 0)
        error('Atom coordinates in the first lattice vector direction are out of cell');
    end
else
    S.Atoms(S.Atoms(:,1) >= S.L1,1) = S.Atoms(S.Atoms(:,1) >= S.L1,1) - S.L1;
    S.Atoms(S.Atoms(:,1) < 0,1) = S.Atoms(S.Atoms(:,1) < 0,1) + S.L1;
    disp('INFO: atom positions are changed in my_check_atomlocation')
end

if S.BCy == 1
    if (sum(S.Atoms(:,2) >= S.L2 | S.Atoms(:,2) < 0) > 0)
        error('Atom coordinates in the second lattice vector direction are out of cell');
    end
else
    S.Atoms(S.Atoms(:,2) >= S.L2,2) = S.Atoms(S.Atoms(:,2) >= S.L2,2) - S.L2;
    S.Atoms(S.Atoms(:,2) < 0,2) = S.Atoms(S.Atoms(:,2) < 0,2) + S.L2;
    disp('INFO: atom positions are changed in my_check_atomlocation')
end

if S.BCz == 1
    if (sum(S.Atoms(:,3) >= S.L3 | S.Atoms(:,3) < 0) > 0)
        error('Atom coordinates in the third lattice vector direction are out of cell');
    end
else
    S.Atoms(S.Atoms(:,3) >= S.L3,3) = S.Atoms(S.Atoms(:,3) >= S.L3,3) - S.L3;
    S.Atoms(S.Atoms(:,3) < 0,3) = S.Atoms(S.Atoms(:,3) < 0,3) + S.L3;
    disp('INFO: atom positions are changed in my_check_atomlocation')
end

end % end function
