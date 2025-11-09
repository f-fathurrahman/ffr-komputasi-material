%copyright by J. E Hasbun and T. Datta
%qholev.m
%Plots the harmonic oscillator potential and a quantum
%mechanical wavefunction for n=2.
%Dimensionsless unit y is used, such that x=sqrt(hbar/(m*w))y,
%Energy is in units of hbar*omega/2. Omega is the natural frequency.
clear;
dy=2*3/100;               %step size
y=-3:dy:3;                %range
H(0+1,:)=ones(1,101);     %Hermite polynomials for n=0,1
H(1+1,:)=2*y;
H(2+1,:)=2*y.*H(1+1,:)-2*1*H(0+1,:); %recursion formula
Expo=exp(-y.^2/2);
V=@(y) y.^2;              %potential function in dimensionless units
plot(y,V(y),'k:')
hold on
eps=2*2+1;
A=1/sqrt(2^2*factorial(2)*sqrt(pi));
psi=A*H(2+1,:).*Expo;
plot(y,psi+eps,'k')       %add energy level to psi
xlabel('y')
ylabel('\psi(y),V(y)')
