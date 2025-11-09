%copyright by J. E Hasbun and T. Datta
%group_phase_vel_example.m
%This code simulates the displacent the wave y1=A*sin(k1*x-w1*t)
%as a function of space and time.
clear;
A=1.0;
N=401;                        %x points
w1=10;                        %frequency
k1=0.1;                       %wavevector
tau1=2*pi/w1;                 %period
tmax=tau1;                    %maximum time
lambda1=2*pi/k1;              %wavelength
lmax=2*lambda1;               %maximum distance
vph1=w1/k1;
fprintf('w1=%5.3f, k1=%5.3f, vph1=w1/k1=%5.3f\n',w1,k1,vph1)
x=0:lmax/N:lmax;              %position array
t=0:0.01:2*tmax;              %time array
for i=1:length(t)             %loop over the time array
  clf                         %clear figure before replotting
  y1=A*sin(k1*x-w1*t(i));     %calculate the wave over all space
  plot(x,y1,'k');             %replot y1 for each time
  axis([0 lmax -A A])
  pause (0.05)                %pause momentarily to see figure
end
xlabel('x (m)')
ylabel('y_1(x,t)')
