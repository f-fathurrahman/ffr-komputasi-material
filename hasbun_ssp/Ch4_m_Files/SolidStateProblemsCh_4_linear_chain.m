%copyright by J. E Hasbun and T. Datta
%linear_chain.m
%simulation of the time dependent motion of N+1 particle
%oscillator chain model. We use copper atoms equidistantly displaced
%from each other by the lattice constant a.
%We take initial distribution as yA*sin[j*pi/L], where L=N*a to be the
%initial distribution of particles. Use M=1.05e-25kg for the mass of a
%Copper atom, and run the simulation for a short time enough to
%determine the oscillation period from which w can be determined.

function linear_chain
clear
global yA L C Np1 M;
a=2.56e-10;       %lattice constant in Copper in m
yA=1.0;           %amplitude in Angstroms
N=25;             %number of cells;
Np1=N+1;          %number of particles
L=N*a;            %length over which particles are spread
B=1.42e11;        %Bulk modulus for Copper in N/m^2
C=B*a;            %spring constant in Copper
if(Np1<3), return; end;    %no less than three particles allowed
x0=0:a:L;         %x positions for N particles, evenly space
M=1.05e-25;       %Copper atom mass in kg
nm=1;             %the nm-th mode
t0=0.0;           %set time range
tf=3.75e-12;      %seconds
range=[t0;tf];
y0=profile(x0,nm);       %Np1 particles y initial position (row vector)
ymax=max(y0);
v0=zeros(1,Np1);         %particle not moving initially (row vector)
ic=[y0';v0'];            %the initial conditions, as columns
%opt=odeset('AbsTol',1.e-7,'RelTol',1.e-7); %if needed use this line
%[t,y]=ode23(@derivs,range,ic,opt);         %with this one, but comment the
                                            %next one
[t,y]=ode23(@derivs,range,ic);
%Plot the height of the particles at each x position every time
for i=1:length(t)
    %fprintf('t=%5.2f\n',t(i));
    %fprintf('%6.2f',y(1:N,i));
    %disp(' ');
    clf
    plot(x0,y(i,1:Np1),'k.','MarkerSize',5)
    axis([x0(1) x0(Np1) -ymax ymax])
    str1=cat(2,'t=',num2str(t(i),'%5.2e sec'));
    text((x0(Np1)-x0(1))*0.01,ymax*(1-0.1),str1,'FontSize',14);
    pause(0.05)
end
xlabel('x (m)','FontSize',14)
ylabel('y (Angstroms)','FontSize',14),
title('Oscillator Chain Model','FontSize',14)

function y=profile(x,nm)
global yA L;
%Possible initial height profile of the particles
y=yA*sin(nm*x*pi/L);    %nth vibrational mode

function [der]=derivs(t,y)
global C Np1 M;
%This function evaluates all the components of the
%derivative vector, putting the results in the array 'der'.
%y(1:Np1): y positions for the N+1 particles
%y(Np1+1:2*Np1): velocities for the N+1 paricles
der=zeros(2*Np1,1);    %initialize der as a column vector
for i=2:Np1-1
    der(i) =y(Np1+i);  %velocities
end
der(1)=0.0;            %end particles are fixed,so
der(Np1)=0.0;          %their y velocities are zero.
for i=Np1+2:2*Np1-1
    der(i)=C*(y(i-Np1-1)+y(i-Np1+1)-2*y(i-Np1))/M; %accelerations
end
der(Np1+1)=0.0;        %end particles don't move, so
der(2*Np1)=0.0;        %their accelerations are zero
