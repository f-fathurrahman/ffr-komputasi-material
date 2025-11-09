%Example code to carry out the conversion process from crystal
%axes to Cartesian axes.
clear                 %It is a good idea to clear the memory
a=                    %axes in angstroms
b=
c=
alpha=                %angles
beta=
gamma=
%Define c1, c2, and c3 to compose the C matrix
c1= c*cos(beta)
c2= c*(cos(alpha)-cos(gamma)*cos(beta))/(sin(gamma))
c3= + sqrt((c^2 - (c1)^2 - (c2)^2))
C =[[a b*cos(gamma) c1] [0 b*sin(gamma) c2] [0 0 c3 ]]  %C matrix
%The coefficients of the atoms in the crystal representation
u=
v=
w=
%Print the vector in the Cartesian representation
C*[u;v;w]
