%
% construct the SiH4 (Silane) molecule 
%
% 1. construct atoms
%
a1 = Atom('Si');
a2 = Atom('H');
atomlist = [a1 a2 a2 a2 a2];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
redxyz = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzlist = redxyz*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut', 5.0, 'name','SiH4' );

fprintf('Done setting up SiH4.\n')