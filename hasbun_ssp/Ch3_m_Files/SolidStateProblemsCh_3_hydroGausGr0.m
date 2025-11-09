%copyright by J. E Hasbun and T. Datta
%hydroGausGr0.m
%This script is set up for two gaussians.
%Program to find the ground state if Hydrogen based on expanding
%the ground state wavefunction with Gaussian orbitals
%(psi(r)=sum_i Ci*Gi(r)). The eigenvalue problem is
%[-(hbar^2/(2*m)) del^2 -ke^2/r]psi=E*psi. If we let r=rbar*a0 and
%E=Ebar*Eb, where a0=4pi*epsilon0*hbar^2/me^2 and
%Eb=2*hbar^2/(2m*a0^2)=1 Hartree (atomic units), the SE becomes
%[-del^2/2-1/r]psi=E*psi with E in units of Eb, and r in units of a0.
%The ground state <H>=<psi|H|psi>/<psi|psi> is calculated as
%<H>=sum(all of Ci*Hij*Cj)/sum (all Ci*Sij*Cj) or just multiply matrices
%once all matrix elements and overlaps have been found.
%Expected result is 1/2 a Hartree.
%Matlab's fminsearch minimizes the energy & optimizes the Gaussian
%exponents, as well the C's to have a total of 2*nG parameters to
%optimize, where nG is the number of Gaussians employed.

function hydroGausGr0
clear all
global Cg
e=1.602176487e-19;     %electronic charge
h=6.62606896e-34;      %Planck'constant (J.s)
eps0=8.854187817e-12;  %Permittivity of free space (C^2/N/m^2)
k=1/4./pi/eps0;        %Electrical constant (N.m^2/C^2)
hbar=h/2./pi;          %hbar
me=9.10938215e-31;     %electron mass (kg)
a0=hbar^2/me/e^2/k;    %Bohr radius (m)
Eb=2*hbar^2/(2*me*a0^2);
%The Gaussian exponents guesses (play with them)
%A way to guess them is to fit the Hydrogen ground state also.
alphaG(1)=1.7792e-001;
alphaG(2)=1.9701e+000;
nG=length(alphaG);
Cg(1:nG)=1/nG;         %starting guesses for the C's as similar weights
%Put the alphas into half of ParsG and the C's in the other half
ParsG=zeros(1,2*nG);   %declare the ParsG array size
ParsG(1:nG)=alphaG(1:nG);
ParsG(nG+1:2*nG)=Cg(1:nG);
%
opts =optimset('TolFun',1e-10,'TolX',1.e-10,'MaxIter',800);
[ParsF,Ene,Eflag,Output]=fminsearch(@EneGroundFind2,ParsG,opts);
disp('Results')
fprintf('Energy=%9.6f (Hartree) or %9.6f (eV)\n',Ene,Ene*Eb/e)
disp('The alphas:')
fprintf('alpha_Guess=['), fprintf('%9.4e, ',alphaG), disp(']') %the guesses
fprintf('alpha_final=['), fprintf('%9.4e, ',ParsF(1:nG)), disp(']')  %the final ones
disp('The C''s:')
fprintf('C_Guess=['), fprintf('%9.4e, ',Cg), disp(']') %the guesses
fprintf('C_final=['), fprintf('%9.4e, ',ParsF(nG+1:end)), disp(']')  %the final ones
fprintf('iterations=%4i\n',Output.iterations)
%disp(Output.message)
%
%Plot the wave function expansion for the Ground state
%The ground state is psi(x)=exp(-x)=sum_over_n (C_n*Gaussian_PSI_n(x))
x=0:0.1:5;
[ana,numG,Error]=WaveG(x,ParsF(1:nG),ParsF(nG+1:end));
fprintf('Estimate of the wavefunction error=%12.6f\n',Error)
plot(x,ana,'k.')
hold on
plot(x,numG,'bo')
hold off
xlabel('r (a_0)')
ylabel('\phi_{1s}(r)')
str=cat(2,num2str(nG),'-Gaussian expansion');
legend('exp(-r)',str)
axis([0 max(x) 0 max(max(numG),max(ana))])

function [ana,numG,Error]=WaveG(x,alpha,Coef)
nG=length(alpha);
ana=exp(-x);
for i=1:length(x)
  numG(i)=0;
  for j=1:nG
    numG(i)=Coef(j)*exp(-alpha(j)*x(i)^2)+numG(i); %Gaussian expansion
  end
end
%Normalization
%Using N^2*4*pi*IntegralOf(r^2 * psi^2 dr), then
%N=1/sqrt(4*pi*IntegralOf(r^2 * psi^2 dr))
dx=(max(x)-min(x))/(length(x)-1);
Cana=1.0/sqrt(dx*trapz(x.^2.*abs(ana).^2))/2/sqrt(pi); %Ground state Norm const
ana=ana*Cana;
CNG=1.0/sqrt(dx*trapz(x.^2.*abs(numG).^2))/2/sqrt(pi); %Gaussian Norm const
numG=numG*CNG;
Error=sqrt(sum(abs((ana-numG).^2)))/length(ana);

function EG=EneGroundFind2(pars)
%1st half of pars are the alphas, 2nd half  are the C's
nG=length(pars)/2;
Over=OverMatrix(pars(1:nG));        %Gaussian Overlap S matrix elements
Dels=DelsMatrix(pars(1:nG));        %-Del^2/2 Gaussian matrix elements
Orints=OrintsMatrix(pars(1:nG));    %1/r Gaussian matrix elements
HM=(Dels+Orints);                   %The total hamiltonian matrix
%easiest to do matrix products: (vector)*Hmatrix*(vector)'
hS=pars(nG+1:2*nG)*HM*pars(nG+1:2*nG)';
OS=pars(nG+1:2*nG)*Over*pars(nG+1:2*nG)';
EG=hS/OS;                           %ground state energy

function SM=OverMatrix(par)
%Gaussians overlap integrals
%See J. M. Thijssen' "computational Physics", chapter 3 (H atom case)
nG=length(par);
for i=1:nG
  for j=i:nG
    SM(i,j)=(pi/(par(i)+par(j)))^1.5;
    SM(j,i)=SM(i,j); %Hermitian
  end
end

function DM=DelsMatrix(par)
%Gaussian integrals for the -Del^2/2 term
%%See J. M. Thijssen' "computational Physics", chapter 3 (H atom case)
nG=length(par);
for i=1:nG
  for j=i:nG
    DM(i,j)=3.*par(i)*par(j)*pi^1.5/(par(i)+par(j))^(5./2.);
    DM(j,i)=DM(i,j); %Hermitian
  end
end

function OM=OrintsMatrix(par)
%Performs the Gaussian integrals for the -1/r term
%%See J. M. Thijssen' "computational Physics", chapter 3 (H atom case)
nG=length(par);
for i=1:nG
  for j=i:nG
    OM(i,j)=-2*pi/(par(i)+par(j));
    OM(j,i)=OM(i,j); %Hermitian
  end
end
