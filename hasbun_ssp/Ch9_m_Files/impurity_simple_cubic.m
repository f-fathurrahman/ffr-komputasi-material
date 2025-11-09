%copyright by J. E Hasbun and T. Datta
%impurity_simple_cubic.m
%Here we obtain f0=(1/N)sum_over_k (1/(E0-Ek))
%from the density of states. That is,
%f0(E)=int(g(E') dE'/(E-E'+im*delta)) then calculate the
%impurity level from 1/Real(f0)-Emu=0

function impurity_simple_cubic
clear;
delta=1.e-3;
im=complex(0.0,1.0);
a=1.0;             %lattice constant
tpa=2*pi/a;
VBZ=tpa^3;         %total SC BZ volume
e2p=3.0;
e1p=-e2p;
ntpmax=251;
esp=(e2p-e1p)/(ntpmax-1);
%Get the simple cubit density of states to be reused
for ntp=1:ntpmax
  e0p(ntp)=e1p+(ntp-1)*esp;                %E prime
  top(ntp)=jelittoScDosAnal(e0p(ntp))/VBZ; %integration numerator
end
%Search for the impurity bound state
%Use Newton-Raphson x(i+1)=x(i)-dE/(F(E+dE)/F(E)-1), where
%F(E)=1/f0(E)-E_mu is the function whose zero we seek
tol=1.e-3;              %convergence tolerance
Ei=3.25; E_0=0;         %impurity and host diagonal energies
Emu=Ei-E_0;             %The perturbation Emu=E_i-E_0
Eguess=1.2*Emu;         %initial guess (for a deep level)
iFg=10*tol;             %convergence check variable
iter=0;
maxiter=30;             %maximum iterations
%Iteration loop
while(abs(iFg) >= tol & iter < maxiter)
  iter=iter+1;                   %iteration counter
  de=max(0.5,abs(Eguess));       %energy variation
  Eguess_de=Eguess+de;           %vary the energy guess
  for ntp=1:ntpmax
    denog(ntp)=Eguess-e0p(ntp)+im*delta;        %1st denominator
    denog_de(ntp)=Eguess_de-e0p(ntp)+im*delta;  %2nd denominator
  end;
  f0g=real(singInt(top,denog,esp));             %integration
  f0g_de=real(singInt(top,denog_de,esp));       %integration
  iFg=1/f0g-Emu;                %The function whose zero is seek
  iFg_de=1/f0g_de;              %the varied function
  Ecorr=de/(iFg_de/iFg-1.0);    %Newton-Raphson correction
  Eguess=Eguess-Ecorr;          %new guess
end
Ebo=Eguess;                     %bound state due to the impurity
disp('Results: note that Ecorr is small at convergence')
fprintf('iter=%4i, Ecorr=%6.4g, 1/f0g=%6.4g, iFg=%6.4g\n',iter,...
  Ecorr,1/f0g,iFg)
fprintf('Peturbation: Emu=%6.4f, Bound state found: Ebo=%6.4f\n',Emu,Ebo)
