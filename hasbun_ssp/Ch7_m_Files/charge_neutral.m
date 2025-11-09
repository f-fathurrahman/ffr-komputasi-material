%copyright by J. E Hasbun and T. Datta
%charge_neutral.m
%It uses the charge neutrality condition to obtain the chemical
%potential (mu) self-consistently based on the donor (ND) and
%acceptor (NA) concentration. Once mu is obtained, the electron
%concentration (n) and the hole (p) concentration can be calculated
%at a given temperature.
clear;
hbarC=197;     %hbar*c in units of eV*nm
kB=8.6174e-5;  %Boltzmann constant in units of eV/K
Ev=0.0;        %valence band edge in eV
Ec=1.0;        %conduction band edge in eV
Eg=Ec-Ev;      %band gap
ED=0.050;      %donor binding energy in eV
EA=0.025;      %acceptor binding energy in eV
epsd=Ec-ED;    %donor impurity level
epsa=Ev+EA;    %acceptor impurity level
mh=0.025;      %hole effective mass value in units of me
mc=0.05;       %electron effective mass value in units of me
meCs=0.5e6;    %electron mass in eV
T=300;         %Kelvin
tau=kB*T;
ND=1e24;       %donor concentration in 1/m^3
NA=1e14;       %acceptor concentration in 1/m^3
%Actual calculations use units of distance in nm
NDnm=ND*(1e-9)^3; %converting 1/m^3 to 1/nm^3
NAnm=NA*(1e-9)^3; %as used here in the FN formula
mu0=Ev+Eg/2+(3/4)*tau*log(mh/mc); %intrinsic mu - use as guess
%Based on the neutrality condition, the function FN must be
%zero for the correct chemical potential mu (variable x).
n0=2*(mc*meCs*tau/hbarC^2/2/pi)^(3/2);
p0=2*(mh*meCs*tau/hbarC^2/2/pi)^(3/2);
FN =@(x) n0*exp((x-Ec)/tau)+NAnm/(1+2*exp((epsa-x)/tau))- ...
    p0*exp((Ev-x)/tau)-NDnm/(1+2*exp((x-epsd)/tau));
%
mu=fzero(FN,mu0); %solve for mu using mu_0 as guess
%use the new mu to get the electron and hole concentrations
nnm=n0*exp((mu-Ec)/tau); %electron conc.
pnm=p0*exp((Ev-mu)/tau); %hole conc.
n=nnm/(1e-9)^3; %converting to 1/m^3
%The approximations for when NA=0
n_D_app=sqrt(n0*NDnm/2)*exp(-ED/tau/2); %elect. conc.
n_D_approx=n_D_app/(1e-9)^3; %converting to 1/m^3
mu_approx=Ec-ED/2+(tau/2)*log(NDnm/2/n0);%chem. pot.
fprintf('mh=%6.3f me, mc=%6.3f me, T=%6.2f K\n',mh,mc,T)
fprintf('Ev=%6.3f eV, Ec=%6.3f eV, mu0=%6.3f eV\n',Ev,Ec,mu0)
fprintf('NA=%6.3g 1/m^3, ND=%6.3g 1/m^3, mu=%6.3f eV\n',NA,ND,mu)
fprintf('n=%6.3g 1/m^3\n',n)
disp('Next are the approximate results assuming NA=0')
fprintf('n_approx=%6.3g 1/m^3, mu_approx=%6.3f eV\n',n_D_approx,mu_approx)
