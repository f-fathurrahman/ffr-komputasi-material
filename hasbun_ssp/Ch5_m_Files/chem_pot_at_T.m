%copyright by J. E Hasbun and T. Datta
%chem_pot_at_T.m
%Chemical potential at a given temperature. The integral
%of the density of orbitals is integrated over energy.
%The Newton-Raphson method is used to get the chemical potential
%self-consistently. The integration is done by the trapezoidal rule.
function chem_pot_at_T
clear
global fe De
h=6.62606896e-34;           %Planck'constant (J.s)
hbar=h/2./pi;               %hbar
me=9.10938215e-31;          %electron mass (kg)
e=1.602176487e-19;          %electronic charge
Eb=e;                       %energy in Joules (=1eV)
kB=1.3806504e-23;           %Boltzmann Constant (J/K)
Tb=Eb/kB;                   %temperature unit
Lb=(hbar^2/2/me/Eb)^(1/2);  %length unit
fprintf('Energy unit: %3.2f eV,\n',Eb/e)
fprintf('Temperature unit: %4.3g K\n',Tb)
fprintf('Length unit: %4.3g m\n',Lb)
Ef=2.5;            %Fermi level at T=0 (in Eb units)
N=Ef^(3/2)/3/pi^2; %Electron number corresponding to T=0
N_over_V=(2*me/hbar^2)^(3/2)*(Ef*Eb)^(3/2)/3/pi^2; %actual N/V
fe= @(eps,mu,T) 1./(exp((eps-mu)./T)+1);  %FD distribution definition
De= @(eps) eps.^(1/2)/2/pi^2;             %Density of states 3D definition
fprintf('T=0 case: Ef=%3.2g, N=%4.3f\n',Ef,N)
N_iter=0;          %iteration counter
N_iterMax=5;       %maximum iteration number
mu_old=Ef;         %initial guess is the Fermi level at T=0
mu_new=mu_old;
del=1.e-3;
mu_corr=10*del;    %use the correction for tolerance
mu_max=8.0*Ef;     %use a large number for the integral's upper limit
Es=(mu_max-0)/500; %even divisions=> odd variables
eps=0:Es:mu_max;
%Let's find how the Fermi level changes with temperature next
%temperature variable in units of Tb
T=0.5;             %temperature variable
while (abs(mu_corr) > del & N_iter < N_iterMax)
  mu1=mu_new;
  mu2=mu1+del;             %vary mu by small ammount
  y1=norb(eps,mu1,T);      %calculate with mu
  y2=norb(eps,mu2,T);      %calculate with mu+del
  %For integration use the trapezoidal
  N1=trapz(eps,y1);        %integral=> occupied states with mu
  N2=trapz(eps,y2);        %integral=> occupied states with mu+del
  ff=N-N1;                 %function whose zero we seek at mu
  ffd=N-N2;                %function at mu+del
  %Newton Raphson: x(i+1)=x(i)-f(xi)/f'(xi), where f'~(f(x+del)-f(x))/del
  %or x(i+1)=x(i)-del/(f(x+del)/f(x)-1)
  N_iter=N_iter+1;
  mu_corr=del/(ffd/ff-1.0);
  mu_new=mu_old-mu_corr;   %Newton-Raphson step
  mu_old=mu_new;
end
fprintf('T (Tb)=%4.2f, N_iter=%3i\n',T,N_iter)
fprintf('ff=%3.2g, mu_corr=%3.2g, result_N=%4.3f\n',ff,mu_corr,N1)
fprintf('N/V (1/m^3)=%5.3e, final mu (eV)=%5.3f\n',N_over_V,mu_new)


function y=norb(eps,mu,T)
global fe De
for i=1:length(eps)
  y(i)=fe(eps(i),mu,T)*De(eps(i));
end
