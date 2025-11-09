%copyright by J. E Hasbun and T. Datta
%thermal_potential.m
%The script plots the Lennard-Jones (LJ) potential and it's expanded
%approximation to 4rth order in (r-R). We use energy/atom in
%units of epsilon, and distance in units of sigma.
%It also calculates the average

function thermal_potential
clear
global alpha U0
alpha=(2*12.13/14.45)^(1/6);  %LJ bond length in Sigma units
U0=14.45^2/2/12.13;           %LJ min energy magnitude in epsilon units
rl=0.1; rs=0.001; ru=15;
r=rl:rs:ru;
plot(r,ULJ(r),'k','LineWidth',2)
hold on
plot(r,ULJexp(r,0),'k:','LineWidth',2)
plot(r,ULJexp(r,1),'k--')
h=legend('Lennard-Jones','4th-Order Approx','Harmonic','Location','North');
set(h,'FontSize',14)
axis ([alpha-0.4 alpha+0.4 -U0*(1+0.4) 5*U0])
xlabel('r (\sigma)','FontSize',14)
ylabel('U(r) (\epsilon)','FontSize',14)

function Ut=ULJ(x)
%Lennard-Jones potential function.
%Energy in units of epsilon, distance in units of sigma.
Ut=2*(12.13./x.^12-14.45./x.^6);

function Ua=ULJexp(x,k)
global alpha U0
%Fourth order Lennard-Jones potential expanded function.
%Energy in units of epsilon, distance in units of sigma.
Cpp=12.13*12*13/alpha^14-14.45*6*7/alpha^8;
Cppp=-12.13*12*13*14/alpha^15+14.45*6*7*8/alpha^9;
Cpppp=12.13*12*13*14*15/alpha^16-14.45*6*7*8*9/alpha^10;
if (k==0)                      %up to 4rth order in x-alpha
  Ua=-U0+2*Cpp*(x-alpha).^2/2+2*Cppp*(x-alpha).^3/6+...
    2*Cpppp*(x-alpha).^4/24;
else
  Ua=-U0+2*Cpp*(x-alpha).^2/2; %the quadratic case only
end
