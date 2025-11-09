%copyright by J. E Hasbun and T. Datta
%hydro_mol_ion_sym.m
%This program incorporates the analytic solution to the
%singly ionized hydrogen molecule.
%
clear
e=1.602176487e-19;       %electronic charge
h=6.62606896e-34;        %Planck'constant (J.s)
eps0=8.854187817e-12;    %Permittivity of free space (C^2/N/m^2)
k=1/4./pi/eps0;          %Electrical constant (N.m^2/C^2)
hbar=h/2./pi;            %hbar
me=9.10938215e-31;       %electron mass (kg)
a0=hbar^2/me/e^2/k;      %Bohr radius (m)
Eb=2*hbar^2/(2*me*a0^2); %Hartree energy unit
%Analytic result
R=0.25:0.001:5;                  %plotting range in unit of a0
Ab=1./R-(1./R+1.0).*exp(-2*R);   %A in units of Eb
Bb=(R+1).*exp(-R);               %B in units of Eb
Del=(1+R+R.^2/3).*exp(-R);       %overlap, R in a0 units
E1s=-1.0/2.0;                    %Bohr energy in Eb units
Eions=1./R;                      %Ion-Ion repulsion energy
Ebo=E1s-(Ab+Bb)./(1.0+Del);      %Bonding energy
Etotbo=Ebo+Eions;                %Total energy for bonding
[Etmin,index_ymin]=min(Etotbo);  %total energy minimum & its array index
Rf=R(index_ymin);                %equilibrium bond length at min of energy
Ebomin=Ebo(index_ymin);          %bonding energy at equilibrium
Eimin=Eions(index_ymin);         %ion energy at equilibrium
%Plots follow
plot(R,Ebo,'k-.')
hold on
plot(R,Eions,'k:')
plot(R,Etotbo,'k-')
line([0 max(R)],[E1s E1s],'Color','k','LineStyle','--') %H atom energy
xlabel('R(a_0)'), ylabel('E(\epsilon_b)')
axis([0 max(R) min(Ebo) -min(Ebo)])
legend('<E>','E_{ions}','E total','atomic: E1s')
title('Ionized hydrogen molecule: bonding state')
hold off
%
disp('H2+ molecule analytic Results - bonding state')
fprintf('R_equilibrium=%9.4f, (a0) or %9.6f (Angs)\n',Rf,Rf*a0/1e-10)
fprintf('<E>=%9.6f (Hartree) or %9.6f (eV)\n',Ebomin,Ebomin*Eb/e)
fprintf('Eions=%9.6f (Hartree) or %9.6f (eV)\n',Eimin,Eimin*Eb/e)
fprintf('Etotal=%9.6f (Hartree) or %9.6f (eV)\n',Etmin,Etmin*Eb/e)
