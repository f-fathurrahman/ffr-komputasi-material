%copyright by J. E Hasbun and T. Datta
%ionic_NaCl.m
%This script performs the calculation of the potential
%U=z*lambda*exp(-R/rho)-UM*R0./R, where UM=alpha*k*e^2/R0 is
%the Madelung energy. The NaCl ionic system is done here.
%It makes use of the Madelung constant for the FCC structure.
%It calculates the minimum energy (lattice energy), U0,
%and the nearest neighbor distance, R0.
clear
Na=6.02214179e23;       %Avogadro's constant (1/mol)
JpK=4.186e3;            %Joules per Kcal
e=1.602176487e-19;      %electronic charge
eps0=8.854187817e-12;   %Permittivity of free space (C^2/N/m^2)
k=1/4/pi/eps0;          %constant (mks units)
alpha=1.747;            %NaCl Madelung Constant parameter
rho=0.321e-10;          %NaCl potential decay parameter (m)
Zlamb=6553.59;          %NaCl z*lambda parameter (ev)
const=k*e^2*alpha/rho/(Zlamb*e);        %dimensionless, to get R0
fR0=inline('x.^2.*exp(-x)-c','x','c');  %function=0 to get R0
xg=2-log(2^2*exp(-2))-log(const);       %rough guess for x
[x,fval] =fzero(fR0,xg,[],const);       %finds x as the zero of fR0
R0=x*rho;               %R0 in meters
%R0=2.82e-10            %table value for R0 if desired (m) for NaCl
UM=alpha*k*e^2/R0/e;    %Madelung energy per molecule in eV
U0=-UM*(1-rho/R0);      %U0  in eV
U0_KpM=U0*e*Na/JpK;     %U0(minimum energy)in Kcal per mol
%Define U_of_R per molecule in eV
U_of_R=@(R,Zlamb,rho,UM,R0) Zlamb*exp(-R/rho)-UM*R0./R;
Rl=0.6*R0; Ru=3*R0;
R=Rl:(Ru-Rl)/50:Ru;           %R range
U=U_of_R(R,Zlamb,rho,UM,R0);  %U(R) per molecule in eV
h=plot(R*1e10,U,'k');
hold on
plot(R0*1e10,U0,'k*')
xlabel('R (Angstroms)')
ylabel('U (eV)')
disp('NaCl ionic system')
fprintf('Madelung energy=%8.4f eV, R0=%8.4f Angstroms\n',UM,R0/1e-10)
fprintf('Lattice energy=%9.4f eV =%10.4f Kcal/mol\n',U0,U0_KpM)
str=cat(2,'NaCl potential: \rho = ',num2str(rho/1e-10,'%6.4f'),...
  ' Ang, z\lambda = ',num2str(Zlamb,'%4.2f'),' eV');
title(str)
str2=cat(2,'R_0 = ',num2str(R0/1e-10,'%7.2f'),...
  ' Ang, U_0 = ',num2str(U0,'%7.2f'),' eV, =',...
  num2str(U0_KpM,'%8.2f'),' Kcal/mol');
legend('U(R) for NaCl',str2,0)
