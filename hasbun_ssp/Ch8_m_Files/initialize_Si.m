%copyright by J. E Hasbun and T. Datta
%initialize_Si.m
function [a,esc,esa,epc,epa,ess,esp,exx,exy,system]=initialize_Si()
system='Si';
%The off diagonal coefficients
esssig=-1.40;
espsig=1.84;
eppsig=3.24;
epppi=-0.81;
h=6.626075e-34;   %Planck's constant in J-sec
hbar=h/2/pi;
me=9.1093897e-31; %electron mass (kg)
e=1.60217733e-19; %electron charge (C). Also recall 1 Joule=(1/e)eV
const=(hbar^2/me)*(1/e)*(1e10)^2; %hbar^2/me in eV-Angtrom^2=7.62eV-Angs^2
%=============== semiconductor parameters ===========================
%Paramaters from Harrison's Electronic Structures and the Properties of solids
%Elements: cation=anion (eV), bond length (angstroms) inputs
%c=cation, a=anion, es=s-energy, ep=p-energy
esa=-13.55; epa=-6.52; esc=esa; epc=epa; d=2.35; %Silicon
a=4*d/sqrt(3);     %Lattice constant
%zinc-blende ion near neighbor positions are located at
%{[1,1,1]a/4, [1,-1,-1]a/4, [-1,1,-1]a/4, [-1,-1,1]a/4, so let scale=a/4
scale=a/4;         %near neighbor position vector magnitudes
ro=scale*sqrt(3.); %internuclear or bond length=sqrt(1^2+1^2+1^2)*a/4
r1=ro^2;
vsssig=const*esssig/r1;
vspsig=const*espsig/r1;
vppsig=const*eppsig/r1;
vpppi=const*epppi/r1;
ess=vsssig;
esp=vspsig/sqrt(3);
exx=vppsig/3.+(2./3.)*vpppi;
exy=(vppsig-vpppi)/3.;
