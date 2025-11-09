%copyright by J. E Hasbun and T. Datta
%diagHamil.m
function diagHamil(esc,esa,epc,epa)
%c=cation, a=anion, es=s-energy, ep=p-energy
%builds the diagonal part of the 8x8 Harrison hamiltonian
global H
H(1,1)=esc;
H(2,2)=esa;
H(3,3)=epc;
H(4,4)=epc;
H(5,5)=epc;
H(6,6)=epa;
H(7,7)=epa;
H(8,8)=epa;
