%copyright by J. E Hasbun and T. Datta
%thermal_latt_const.m
%This script calculates the average lattice constant for Argon
%and compares with the experimental values.

function thermal_latt_const
clear
global alpha U0       %variables to be recognized by functions
e=1.602176487e-19;    %electronic charge
kB=1.3806504e-23;     %Boltzmann Contant  (J/K)
%Data for solid Argon digitized from
%"Measurements of X-Ray Lattice Constant,
%Thermal Expansivity, and Isothermal Compressibility
%of Argon Crystals"  O. G. Peterson, T D. N. Batchelder,
%& R. O. Simmons, Phys. Rev. vol. 150, No. 2 (1966)
%T_e is in Kelvin, a_e is Argon's lattice constant converted to angs.
T_e=[4.516,14.837,19.802,20.881,24.762,27.999,29.727,...
  34.910,39.665,49.405,49.627,59.586,60.446,61.750,...
  62.631,63.476,66.961,69.123,70.862,74.976,79.097,...
  79.966,83.006];
a_e=[5.311,5.310,5.317,5.318,5.322,5.327,5.329,5.337,...
  5.346,5.370,5.372,5.397,5.397,5.402,5.410,5.404,...
  5.420,5.425,5.432,5.442,5.455,5.458,5.469];
alpha=(2*12.13/14.45)^(1/6); %LJ bond length in Sigma units
U0=14.45^2/2/12.13;          %LJ min energy magnitude in epsilon units
%Let the extrapolated value at T=0 be the a0 for Argon, then
%R0=a0/sqrt(2) and since R0/sigma=alpha, we find that sigma=R0/alpha
p=polyfit(T_e(1:8),a_e(1:8),1); % fit low T data with a quadratic
a0=polyval(p,0);  %extrapolate to get the T=0 lattice constant;
R0=a0/sqrt(2);    %T=0 bond length for Argon (fcc nearest neighbor)
sigma=R0/alpha;   %our sigma - Lennard-Jones distance parameter
TD=92;            %argon Debye temperature
eps=kB*TD/e;      %energy unit (eV) - Lennard-Jones energy parameter
fprintf('[a0,R0]=[%6.3f, %6.3f]Angs, sigma=%6.3f Angs\n',a0,R0,sigma)
fprintf('Debye-T=%6.2f K, epsilon(eV)=%6.3g eV\n',TD,eps)
%Integrals are evaluated for energy units of epsilon, distance in
%units of sigma and temperature in units of argon's Debye temperatue.
%After integration the conversion is done for temperature and distance
rl=0.1; rs=0.001; ru=15; %integration range in units of sigma
r=rl:rs:ru;
T=0.02:0.01:0.9;
for i=1:length(T)
  %integral based on the fourth order expansion
  yLJexp1=rs*trapz(r,r.*Boltz(ULJexp(r),T(i)));
  yLJexp2=rs*trapz(r,Boltz(ULJexp(r),T(i)));
  raveLJexp(i)=yLJexp1/yLJexp2;     %4rth order approx average bond length
end
%Simple analytic approx
Cpp=12.13*12*13/alpha^14-14.45*6*7/alpha^8;
Cppp=-12.13*12*13*14/alpha^15+14.45*6*7*8/alpha^9;
const=abs(Cppp/Cpp^2)/4;
%for the simple approx add alpha: the equilibrium point in sigma units
r_simple=const*T+alpha;
%For the FCC: lattice constant=sqrt(nearest neighbor distance)
raveLJexp=raveLJexp*sqrt(2);
r_simple=r_simple*sqrt(2);
%Convert temperature to Kelvin and distance to angstroms
T=T*TD;
raveLJexp=raveLJexp*sigma;
r_simple=r_simple*sigma;
hold on
plot(T,raveLJexp,'k--') %Full fourth order expansion case
plot(T,r_simple,'k:')   %Simple analytic result
plot(T_e,a_e,'ko')      %Ar data
legend('<r> (4rth order)','<r> (simple approx)','Ar Data',2)
axis ([0 T(end) a_e(1)*(1-0.005) a_e(end)*(1+0.005)])
xlabel('T (Kelvin)','FontSize',14)
ylabel('Lattice Constant (Angstroms)','FontSize',14)
str=cat(2,'Argon: \epsilon=',num2str(eps,'%6.3g'),...
  'eV, \sigma=',num2str(sigma,'%6.3f'),' Angstroms');
title(str,'FontSize',14);


function y=Boltz(x,T)
%The Boltzmann factor
y=exp(-x/T);

function Ua=ULJexp(x)
global alpha U0
%Fourth order Lennard-Jones potential expanded function.
%Energy in units of epsilon, distance in units of sigma.
Cpp=12.13*12*13/alpha^14-14.45*6*7/alpha^8;
Cppp=-12.13*12*13*14/alpha^15+14.45*6*7*8/alpha^9;
Cpppp=12.13*12*13*14*15/alpha^16-14.45*6*7*8*9/alpha^10;
%up to 4rth order in x-alpha
Ua=-U0+2*Cpp*(x-alpha).^2/2+2*Cppp*(x-alpha).^3/6+...
  2*Cpppp*(x-alpha).^4/24;
