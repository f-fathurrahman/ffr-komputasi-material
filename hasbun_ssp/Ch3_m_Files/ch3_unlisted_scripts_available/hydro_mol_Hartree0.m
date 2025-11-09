%hydro_mol_Hartree0.m by J. E. Hasbun (5/2011)
%This script is set up for two gaussians.
%This program incorporates the wavefunction Gaussian expansion.
%We vary the distance between the two protons to get the potential
%and the energies at the minimum of the potential.
%phi0,1(r)=(s1(r-Ra) +(-) s2(r-Rb))/sqrt(2). Where Ra and Rb are the 
%positions of the two protons. 
%The Hartree approximation makes the ansatz that for the 
%ground state Psi_H0=phi0(r1)*phi0(r2), and
%for the antibonding state Psi_H1=phi1(r1)*phi1(r2).
%The functions used here are listed as follows:
%(a) overlap: Overlap1s(alpa,Ra,beta,Rb);
%(b) Kinetic: Dels1s(alpa,Ra,beta,Rb);
%(c) Coulomb attract proton-electron: Orint1s(alpa,Ra,beta,Rb,Z,Rc);
%(d) Coulomb repel electron-electron: EeInt1s(alpa,Ra,beta,Rb,gama,Rc,delt,Rd).

function hydro_mol_Hartree0
clear all
global isym_ Z Ra 
global ic1_a ic1_b ialpha_a ialpha_b
e=1.602176487e-19;     %electronic charge
h=6.62606896e-34;      %Planck'constant (J.s)
eps0=8.854187817e-12;  %Permittivity of free space (C^2/N/m^2)
k=1/4./pi/eps0;        %Electrical constant (N.m^2/C^2)
hbar=h/2./pi;          %hbar
me=9.10938215e-31;     %electron mass (kg)
a0=hbar^2/me/e^2/k;    %Bohr radius (m)
Eb=2*hbar^2/(2*me*a0^2);
%The Gaussian exponents guesses (obtained by a fit - see exp_Gauss_fit_a,b.m)
%First ion wave function uses alpha
alpha(1)=1.7792e-001;
alpha(2)=1.9701e+000;
nGa=length(alpha);
C1(1:nGa)=1/nGa;       %Real coefficients of expansion all the same

%Second ion wave function uses beta
beta=alpha;
C2=C1;                 %Real coefficients of expansion all the same
%Hydrogen molecule energy with wavefunction of symmetry isym
%where psi(r)=(s1(r-Ra)+isym*s2(r-Rb))/sqrt(2) and si(r-Ri) is the
%ground state wave function expanded in terms of Gaussian exponentials
%H molecule specifics
Z=1.0;               %ion proton number
Ra=0.0;              %Ra - first ion at the origin
Rb_array=0.2:0.01:5; %The variable proton distance
isym=[1,-1]; %symmetries to work with. We do both symmetries below.
% Plot with initial parameters
figure('Name','H2 results with initial parameters')
[EminG,RbG]=plotter2(C1,alpha,C2,beta,isym,Z,Ra,Rb_array);
%Next, we find the minimum energy, bond length, alphas, and betas of 
%the hydrogen molecule numerically.
%We put all the parameters to be varied into array ParsG
%Parameters to be varied, C1's, alpha's, and Rb
disp('Full H2 molecule before minimization');
fprintf('Z=%4.2f, Ra=%9.6f, Rb=%9.6f (a0)\n',Z,Ra,RbG)
fprintf('Etotal=%9.6f (Hartree) or %9.6f (eV)\n',EminG,EminG*Eb/e)
fprintf('alpha_Guess=['), fprintf('%9.4e, ',alpha), disp(']')
fprintf('beta_Guess=['), fprintf('%9.4e, ',beta), disp(']')
fprintf('C1_Guess=['), fprintf('%9.4e, ',C1), disp(']')
fprintf('C2_Guess=['), fprintf('%9.4e, ',C2), disp(']')
%We need to define the parameters to be varied in array ParsG
ic1_a=1; ic1_b=length(C1);
ialpha_a=ic1_b+1; ialpha_b=ic1_b+length(alpha);
ParsG(ic1_a:ic1_b)=C1;                          %C1's to be varied
ParsG(ialpha_a:ialpha_b)=alpha;                 %alpha's to be varied
ParsG(ialpha_b+1)=RbG;                          %Rb to be varied
isym_=isym(1);
pause(1)
disp(' ')
disp('Full H2 molecule after minimization');
opts =optimset('TolFun',1e-6,'TolX',1.e-6,'MaxIter',1600);
[ParsF,Ene,Eflag,Output]=fminsearch(@EFullHmolMin,ParsG,opts);
%translate back to the corresponding parameters being varied
C1=ParsF(ic1_a:ic1_b);
alpha=ParsF(ialpha_a:ialpha_b);
C2=C1;
beta=alpha;
Rb=ParsF(ialpha_b+1);
fprintf('Results for symmetry=%2i\n',isym_);
fprintf('Z=%4.2f, Ra=%9.6f, Rb=%9.6f (a0)\n',Z,Ra,Rb)
fprintf('Etotal=%9.6f (Hartree) or %9.6f (eV)\n',Ene,Ene*Eb/e)
fprintf('alpha=['), fprintf('%9.4e, ',alpha), disp(']')
fprintf('beta=['), fprintf('%9.4e, ',beta), disp(']')
fprintf('C1=['), fprintf('%9.4e, ',C1), disp(']')
fprintf('C2=['), fprintf('%9.4e, ',C2), disp(']')
fprintf('iterations=%4i\n',Output.iterations)
disp(' ')
% Plot with final parameters
figure('Name','H2 results - with final parameters')
[EminF,RbF]=plotter2(C1,alpha,C2,beta,isym,Z,Ra,Rb_array);
%disp('Reading off the final graph one gets these results')
%fprintf('Z=%4.2f, Ra=%9.6f, Rb=%9.6f (a0)\n',Z,Ra,RbF)
%fprintf('Etotal=%9.6f (Hartree) or %9.6f (eV)\n',EminF,EminF*Eb/e)
%
%Summary of results
disp(' ')
disp('Results Summary')
fprintf('Energy=%9.6f (Hartree) or %9.6f (eV)\n',Ene,Ene*Eb/e)
fprintf('R_equilibrium=%9.4f (a0) or %9.6f (Angs)\n',Rb,Rb*a0/1e-10)
for i=1:length(isym)
  [EH21,EH22,Eetot,Etot,Eions]=FullHmol(C1,alpha,C2,beta,isym(i),Z,Ra,Rb);
  if (i==1),
    Ets=Etot;
    disp('Full H2 molecule bonding state: at equilibrium R we have')
  else
    Eta=Etot;
    disp('Full H2 molecule anti-bonding state: at equilibrium R we have')
  end
  fprintf('Z=%4.2f, Ra=%9.6f, Rb=%9.6f (a0)\n',Z,Ra,Rb)
  fprintf('<E1>=%9.6f (Hartree) or %9.6f (eV)\n',EH21,EH21*Eb/e)
  fprintf('<E2>=%9.6f (Hartree) or %9.6f (eV)\n',EH22,EH22*Eb/e)
  fprintf('2<E1>+<E2>=%9.6f (Hartree) or %9.6f (eV)\n',Eetot,Eetot*Eb/e)
  fprintf('Eions=%9.6f (Hartree) or %9.6f (eV)\n',Eions,Eions*Eb/e)
  fprintf('Etotal=%9.6f (Hartree) or %9.6f (eV)\n',Etot,Etot*Eb/e)
  disp(' ')
end
%Final wave functions
figure('Name','H2 wavefunctions with final parameters')
plotterw2(C1,alpha,C2,beta,isym,Z,Ra,Rb,Ets,Eta)

function [EminG,RbG]=plotter2(C1,alpha,C2,beta,isym,Z,Ra,Rb)
for i=1:length(isym)
  for j=1:length(Rb)
    [EH21(j),EH22(j),Eetot(j),Etot(j),Eions(j)]=...
                  FullHmol(C1,alpha,C2,beta,isym(i),Z,Ra,Rb(j));
  end
  if(isym(i)==1)
    [EminG,IG]=min(Etot);  %minimum energy value and its index
    RbG=Rb(IG);            %the value of Rb at min energy (used later)
    %Bonding state
    subplot(1,2,1), plot(Rb,EH21,'k--')
    hold on
    plot(Rb,EH22,'k:')
    plot(Rb,Eetot,'k-')
    plot(Rb,Eions,'k-.')    
    plot(Rb,Etot,'k.','MarkerSize',5)
    line([0 max(Rb)],[-0.5 -0.5],'Color','m','LineStyle','--')
    xlabel('R(a_0)'), ylabel('E(\epsilon_b)')
    Eminsym1=min(Eetot);
    axis([0 max(Rb) Eminsym1 -Eminsym1])
    hl=legend('<E1>','<E2>','2<E1>+<E2>','Eions','E total','atomic: E1s');
    set(hl,'FontSize',8)
    title('Full hydrogen molecule: bonding state')
    hold off
  else
    %Antibonding state
    subplot(1,2,2), plot(Rb,EH21,'k--')
    hold on
    plot(Rb,EH22,'k:')
    plot(Rb,Eetot,'k-')
    plot(Rb,Eions,'k-.')
    plot(Rb,Etot,'k.','MarkerSize',5)
    line([0 max(Rb)],[-0.5 -0.5],'Color','m','LineStyle','--')
    xlabel('R (a_0)'), ylabel('E (\epsilon_b)')
    axis([0 max(Rb) Eminsym1 -Eminsym1])
    hl=legend('<E1>','<E2>','2<E1>+<E2>','Eions','E total','atomic: E1s');
    set(hl,'FontSize',8)
    title('anti-bonding state')
    hold off
  end
end

function plotterw2(C1,alpha,C2,beta,isym,Z,Ra,Rb,Ets,Eta)
%Plots the final wavefunctions using the Gaussian expansion
ru=Rb+1.5*Rb; rl=Ra-1.5*Rb; rs=(ru-rl)/150;
r=rl:rs:ru;
R1=Ra*ones(size(alpha)); %Vector for the first ion position
R2=Rb*ones(size(beta));  %Vector for the 2nd ion position
sq2=1.0/sqrt(2.0);
for i=1:length(isym)
%one way
% phi=C1(1:4)*exp(-(alpha(1:4))'.*((r-Ra).*(r-Ra)))+...
%   isym(i)*C2(1:4)*exp(-(beta(1:4))'.*((r-Rb).*(r-Rb)));
%another way
  for j=1:length(r)
    phi(j)=0.0;
    %Build the Gaussian expansion of the wave functions
    for k=1:length(alpha)
      phi(j)=phi(j)+C1(k)*exp(-alpha(k)*abs(r(j)-R1(k))^2)+...
        isym(i)*(C2(k)*exp(-beta(k)*abs(r(j)-R2(k))^2));
    end
    phi(j)=sq2*phi(j);
  end
  if(isym(i)==1)
    %Bonding state
    subplot(2,1,1), plot(r,phi,'k')
    xlabel('r (a_0)'), ylabel('\phi (2/a_0^{1/3})')
    str=cat(2,'Ra=',num2str(Ra,'%5.3f'),', Rb=',num2str(Rb,'%5.3f'),...
      ' (a0), Etotal=',num2str(Ets,'%6.4f'),' (Hartree)');
    axis([rl ru min(phi)*(1+0.3) max(phi)*(1+0.3)])
    text(rl*(1-0.1),max(phi)*(1+0.1),str,'FontSize',8)
    title('Hydrogen molecule: bonding state')
  else
    %Antibonding state
    subplot(2,1,2), plot(r,phi,'k')
    xlabel('r (a_0)'), ylabel('\phi (2/a_0^{1/3})')
    str=cat(2,'Ra=',num2str(Ra,'%5.3f'),', Rb=',num2str(Rb,'%5.3f'),...
      ' (a0), Etotal=',num2str(Eta,'%6.4f'),' (Hartree)');
    axis([rl ru min(phi)*(1+0.3) max(phi)*(1+0.3)])
    text(rl*(1-0.1),max(phi)*(1+0.1),str,'FontSize',8)
    title('Anti-bonding state')
  end
end

function Etot=EFullHmolMin(pars)
%The first half of the parameters are the C's
%The second half of the parameters are the alphas
global isym_ Z Ra 
global ic1_a ic1_b ialpha_a ialpha_b
%translate the pars to the corresponding parameters being varied
C1=pars(ic1_a:ic1_b);
alpha=pars(ialpha_a:ialpha_b);
C2=C1;
beta=alpha;
Rb=pars(ialpha_b+1);
%Note: the total energy being minimized is Etot=2.0*EH21+EH22+Eions
%or Etot=Eetot+Eions
[EH21,EH22,Eetot,Etot,Eions]=FullHmol(C1,alpha,C2,beta,isym_,Z,Ra,Rb);

function [EneHmol1,EneHmol2,Eetot,Etot,Eions]=...
          FullHmol(C1,alpha,C2,beta,isym_,Z,Ra,Rb)
%The Hartree approximation for the full hydrogen molecule is such that 
%phi0(r)=(s1(r-Ra) + s2(r-Rb))/sqrt(2), 
%phi1(r)=(s1(r-Ra) - s2(r-Rb))/sqrt(2),
%Psi_H0=phi0(r1)*phi0(r2) -> ground state
%Psi_H1=phi1(r1)*phi1(r2) -> anti-bonding state
%This results into two main terms for the energy of the full hydrogen
%molecule. For example, for the ground state we have
%<H(r1,r2)>=2<phi0(r)|Hbo(r)|phi0(r)>/<phi0(r)|phi0(r)> +
%<phi0(r1)*phi0(r2)|h(r1,r2)|phi0(r1)*phi0(r2)>/<phi0(r)|phi0(r)>^2
%Here Hbo is the single electron ionized molecule term and 
%h(r1,r2) is the electron-electron coulom term
%The factor of 2 in the first terms and this two electron term 
%differentiates between the ionized hydrogen molecule and the
%full tw0-electron H molecule. We do the full molecule here.
%We recall the our Si(r-Ri) funtions are expanded in terms of
%Gaussians. That's why we end up summing matrix elements in our approach.
%
R1=Ra*ones(size(alpha)); %Vector for the first ion position
R2=Rb*ones(size(beta));  %Vector for the 2nd ion position
%Overlap S matrix elements, with the symmetry flag included
%We will construct the wavefunction phi(r)=(s1(r)+isym*s2(r))/sqrt(2);
%and find the overlap result
Over=OverMatrixHmol(C1,alpha,C2,beta,R1,R2,isym_);  %The overlaps
Dels=DelsMatrixHmol(C1,alpha,C2,beta,R1,R2,isym_);  %-Del^2/2 Gaussian matrix elements
%1/(r-Ra) Gaussian matrix elements
Orints1=OrintsMatrixHmol(C1,alpha,C2,beta,R1,R2,Z,Ra,isym_);
%1/(r-Rb) Gaussian matrix elements
Orints2=OrintsMatrixHmol(C1,alpha,C2,beta,R1,R2,Z,Rb,isym_);
HM=(Dels+Orints1+Orints2);   %The total hamiltonian matrix
hS=sum(sum(HM));             %Sum all the hamiltonian matrix elements (2D)
OS=sum(sum(Over));           %Sum all the overlap matrix elements (2D)
EneHmol1=hS/OS;              %elect. Ground state energy for the H+ molecule
Eeint=EeIntMatrixHmol(C1,alpha,C2,beta,R1,R2,isym_); %the 4D e-e matrix
%Use function EeIntMatrixHmol2 to check function EeIntMatrixHmol, which
%uses symmetry, and that it works correctly.
%Eeint=EeIntMatrixHmol2(C1,alpha,C2,beta,R1,R2,isym_); %the 4D e-e matrix
%
%Need to sum over all matrix elements to get the whole e-e energy
%contribution - this is a 4-dimensional matrix
Ees=sum(sum(sum(sum(Eeint))));
EneHmol2=Ees/OS^2;                %e-e energy term for the full H2 molecule
Eions=1.0/norm(Rb-Ra);            %Ion-Ion repulsion energy
Eetot=2.0*EneHmol1+EneHmol2;      %total electronic energy
Etot=Eetot+Eions;                 %full H2 molecule total energy

function SM=OverMatrixHmol(C1,alpha,C2,beta,R1,R2,isym_)
%Performs the Gaussians overlap integrals 
%as in J. M. Thijssen text chapter 3 (Hydrogen atom case)
%Construct the wavefunction phi(r)=(s1(r)+isym*s2(r))/sqrt(2)
sq2=1.0/sqrt(2.0);
%Note the coefficients are all real, we don't worry about complex
%conjugates
pars=[alpha,beta];      %exponents vector for respective wavefunction
parC=sq2*[C1,isym_*C2]; %coefficients vector for respective wavefunction
parR=[R1,R2];           %ion positions vector for respective wavefunction
nG2=length(pars);
%Perform the overlap integrals
for i=1:nG2
  for j=i:nG2
    SM(i,j)=Overlap1s(pars(i),parR(i),pars(j),parR(j))*parC(i)*parC(j);
    SM(j,i)=SM(i,j); %Hermitian
  end
end

function DM=DelsMatrixHmol(C1,alpha,C2,beta,R1,R2,isym_)
%Performs the Gaussian integrals for the -Del^2/2 term
%as in J. M. Thijssen text chapter 3 (Hydrogen atom case)
sq2=1.0/sqrt(2.0);
%Note the coefficients are all real, we don't worry about complex
%conjugates
pars=[alpha,beta];      %exponents vector for respective wavefunction
parC=sq2*[C1,isym_*C2]; %coefficients vector for respective wavefunction
parR=[R1,R2];           %ion positions vector for respective wavefunction
nG2=length(pars);
for i=1:nG2
  for j=i:nG2
    DM(i,j)=Dels1s(pars(i),parR(i),pars(j),parR(j))*parC(i)*parC(j);
    DM(j,i)=DM(i,j);    %Hermitian
  end
end

function OM=OrintsMatrixHmol(C1,alpha,C2,beta,R1,R2,Z,Rc,isym_)
%Performs the Gaussian integrals for the -1/r term
%as in J. M. Thijssen text chapter 3 (Hydrogen atom case)
sq2=1.0/sqrt(2.0);
%Note the coefficients are all real, we don't worry about complex
%conjugates
pars=[alpha,beta];      %exponents vector for respective wavefunction
parC=sq2*[C1,isym_*C2]; %coefficients vector for respective wavefunction
parR=[R1,R2];           %ion positions vector for respective wavefunction
nG2=length(pars);
for i=1:nG2
  for j=i:nG2
    OM(i,j)=Orint1s(pars(i),parR(i),pars(j),parR(j),Z,Rc)*parC(i)*parC(j);
    OM(j,i)=OM(i,j); %Hermitian
  end
end

function g=EeIntMatrixHmol(C1,alpha,C2,beta,R1,R2,isym_)
%This function takes advantage if symmetry to perform the 
%Gaussian integrals for the 1/norm(r-r') term
%as in J. M. Thijssen text chapter 3 (Hydrogen atom case)
%The function EeIntMatrixHmol2 does the same but without 
%using symmetry and works the same as this one. Obviously this one
%is preferable for speed purposes.
sq2=1.0/sqrt(2.0);
%Note the coefficients are all real, we don't worry about complex
%conjugates
pars=[alpha,beta];      %exponents vector for respective wavefunction
parC=sq2*[C1,isym_*C2]; %coefficients vector for respective wavefunction
parR=[R1,R2];           %ion positions vector for respective wavefunction
nG2=length(pars);
%The two-electron gaussian matrix elements are calculated
%according to the symmetry rules explained in J. M. Thijssen 
%text chapter 4, pp 72. Here the coefficients parC and the integrals
%are all real.
for ip=1:nG2
  for iq=1:ip     %In these loops, we use symmetry not to run over all
    for ir=1:ip-1 %the basis functions
      for is=1:ir
        %Notice how the iq and ir index labels for g are defined in reverse
        %order to the function call for corresponding parameters
        g(ip,iq,ir,is)=EeInt1s(pars(ip),parR(ip),pars(ir),parR(ir),...
                           pars(iq),parR(iq),pars(is),parR(is))*...
                           parC(ip)*parC(ir)*parC(iq)*parC(is);
        %The following equalities follow symmetry rules (Thijssen)
        g(ir,is,ip,iq)=g(ip,iq,ir,is);
        g(iq,ip,ir,is)=g(ip,iq,ir,is);
        g(ip,iq,is,ir)=g(ip,iq,ir,is);
        g(iq,ip,is,ir)=g(ip,iq,ir,is);
        g(ir,is,iq,ip)=g(ip,iq,ir,is);
        g(is,ir,ip,iq)=g(ip,iq,ir,is);
        g(is,ir,iq,ip)=g(ip,iq,ir,is);
      end
    end
    ir=ip;
    for is=1:iq
      %Again the iq and ir index labels for g are defined in reverse
      %order to the function call for corresponding parameters
      g(ip,iq,ir,is)=EeInt1s(pars(ip),parR(ip),pars(ir),parR(ir),...
                         pars(iq),parR(iq),pars(is),parR(is))*...
                         parC(ip)*parC(ir)*parC(iq)*parC(is);
      %Again, we follow symmetry rules as above (Thijssen)
      g(ir,is,ip,iq)=g(ip,iq,ir,is);
      g(iq,ip,ir,is)=g(ip,iq,ir,is);
      g(ip,iq,is,ir)=g(ip,iq,ir,is);
      g(iq,ip,is,ir)=g(ip,iq,ir,is);
      g(ir,is,iq,ip)=g(ip,iq,ir,is);
      g(is,ir,ip,iq)=g(ip,iq,ir,is);
      g(is,ir,iq,ip)=g(ip,iq,ir,is);
    end
  end
end

function g=EeIntMatrixHmol2(C1,alpha,C2,beta,R1,R2,isym_)
%This function is the equivalent of EeIntMatrixHmol2 but without
%taking advantage of symmetry. It is obviously slower.
%Performs the Gaussian integrals for the 1/norm(r-r') term
%as in J. M. Thijssen text chapter 3 (Hydrogen atom case)
sq2=1.0/sqrt(2.0);
%Note the coefficients are all real, we don't worry about complex
%conjugates
pars=[alpha,beta];     %exponents vector for respective wavefunction
parC=sq2*[C1,isym_*C2]; %coefficients vector for respective wavefunction
parR=[R1,R2];          %ion positions vector for respective wavefunction
nG2=length(pars);
%The two-electron gaussian matrix elements are calculated
%according to the symmetry rules explained in J. M. Thijssen
%text chapter 4, pp 72. Here the coefficients parC and the integrals
%are all real.
for ip=1:nG2
  for iq=1:nG2     %In these loops, we use symmetry not to run over all
    for ir=1:nG2 %the basis functions
      for is=1:nG2
        %Notice how the iq and ir index labels for g are defined in reverse
        %order to the function call for corresponding parameters
        g(ip,iq,ir,is)=EeInt1s(pars(ip),parR(ip),pars(ir),parR(ir),...
          pars(iq),parR(iq),pars(is),parR(is))*...
          parC(ip)*parC(ir)*parC(iq)*parC(is);
      end
    end
  end
end

function y=Overlap1s(alpa,Ra,beta,Rb)
%Performs the Gaussians overlap integrals 
%as in J. M. Thijssen text chapter 3 (1s Gaussian orbital expansion)
%alpa, beta=the alpha and beta are exponents for the Gaussians functions.
%Ra, Rb the atomic nuclei vector positions and could be vectors.
%Note norm(Vector)=sqrt(x^2+y^2+z^2);
ab=(alpa*beta);
apb=alpa+beta;
K=exp(-ab*(norm(Ra-Rb)).^2/apb);
y=K*(pi/apb)^1.5;

function y=Dels1s(alpa,Ra,beta,Rb)
%Performs the Gaussian integrals for the -Del^2/2 term
%as in J. M. Thijssen text chapter 3 (1s Gaussian orbital expansion)
%alpa, beta=the alpha and beta are exponents for the Gaussians functions.
%Ra, Rb the atomic nuclei vector positions and could be vectors.
%Note norm(Vector)=sqrt(x^2+y^2+z^2);
ab=(alpa*beta);
apb=alpa+beta;
abpr=ab/apb;
arg=abpr*(norm(Ra-Rb)).^2;
K=exp(-arg);
y=abpr*(6.0-4.0*arg)*K*(pi/apb)^1.5;
y=y/2;                %the final factor of 1/2 in -Del^2/2;

function y=Orint1s(alpa,Ra,beta,Rb,Z,Rc)
%Performs the Gaussian integrals for the -Z/norm(r-Rc) term
%as in J. M. Thijssen text chapter 3 (1s Gaussian orbital expansion)
%alpa, beta=the alpha and beta are exponents for the Gaussians functions.
%Ra, Rb the atomic nuclei vector positions and could be vectors.
%Note norm(Vector)=sqrt(x^2+y^2+z^2);
ab=alpa*beta;
apb=alpa+beta;
K=exp(-ab*(norm(Ra-Rb)).^2/apb);
Rp=(alpa*Ra+beta*Rb)/apb;
arg=apb*(norm(Rp-Rc))^2;
y=-2*pi*Z*K*F0rpc(arg)/apb;

function y=EeInt1s(alpa,Ra,beta,Rb,gama,Rc,delt,Rd)
%Performs the Gaussian integrals for the 1/norm(r-r') term
%as in J. M. Thijssen text chapter 3 (1s Gaussian orbital expansion)
%alpa, beta, gama, delta=alpha, beta, gamma, delta exponential coefficient
%of the Gaussian functions.
%Ra,Rb,Rc,Rd are the atomic nuclei vector positions and could be vectors.
%Note norm(Vector)=sqrt(x^2+y^2+z^2);
apg=alpa+gama;
bpd=beta+delt;
abgd=1.0/(apg+bpd);
Rp=(alpa*Ra+gama*Rc)/apg;
Rq=(beta*Rb+delt*Rd)/bpd;
K=exp(-alpa*gama*(norm(Ra-Rc)).^2/apg-beta*delt*(norm(Rb-Rd)).^2/bpd);
arg=apg*bpd*abgd*(norm(Rp-Rq)).^2;
y=2.0*pi^2.5*K*F0rpc(arg)*sqrt(abgd)/apg/bpd;

function y=F0rpc(x)
if (abs(x)< 1.e-10)
  y=1.0;
else
  y=sqrt(pi)*erf(sqrt(x))/sqrt(x)/2.0;
end