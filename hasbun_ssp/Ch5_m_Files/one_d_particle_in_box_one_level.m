%copyright by J. E Hasbun and T. Datta
%one_d_particle_in_box_one_level.m
%One-dimensional particle in a box wavefunction solution
%and energy level.
clear
L=0.5;                 %well width in nm
h=6.62606896e-34;      %Planck'constant (J.s)
hbar=h/2./pi;          %hbar
e=1.602176487e-19;     %electronic charge
hbar_eV=hbar/e;        %hbar in eV.s
c=299792458;           %speed of light (m/s)
cnm=c*1e9;             %speed of light (nm/s)
me=0.511003e6/cnm^2;   %electron mass eV.s^2/nm^2
N=301;                 %number of x points to plot
dx=L/(N-1);            %step size for the x variable
%Wave function and energy level follow
hold on                %overlays plots
line([0 0],[0 10],'Color','k')        %left potential wall
line([L L],[0 10],'Color','k')        %right potential wall
E1=pi^2*hbar_eV^2/(2*me*L^2);         %1st energy level in eV
x=0:dx:L;                             %x as an array
psi=sqrt(2/L)*sin(pi*x/L);            %wavefunction versus x array
plot(x,psi+E1,'k')                    %wavefunction plot displaced by En
axis([-0.1 L+0.1 0 10])               %window plot view
xlabel('x (nm)')
ylabel('V (eV),  \psi(x)')
title('One-Dimensional Particle in a Box')
