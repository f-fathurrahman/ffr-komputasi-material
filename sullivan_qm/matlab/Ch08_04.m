% SpecV.m.
% This program calculates the electron density
% in an inﬁnite well with a potential by
% the eigenfunctions and by the spectral method.

clear all
hbar = 1.06e-34;
m0 = 0.25*9.1e-31;
ecoul = 1.6e-19;
eV2J = 1.6e-19;
J2eV = 1./eV2J;
mu = 0.25;
kT = 0.025;

del_x = 2.e-10;

NN = 50;
N_center = NN/2;
% Energies are meV
chi0 = J2eV*hbar^2/(2*m0*del_x^2);

% Convert to nm for plotting
DX = del_x*1e9;
XX = (DX:DX:NN*DX);
mid = DX*N_center;

% This adds the Potential
V = zeros(1,NN);
for n=1:NN
    %V(n) = .004*((N_center-n));
    V(n) = .002*(-n);
end

subplot(3,2,1)
plot(XX,V,'k')
set(gca,'fontsize',12)
grid on
ylabel('V (eV)')
title('SpecV')

% --------- Construct the Hamiltonian ------------
H = zeros(NN,NN);
H(1,1) = 2*chi0+V(1);
H(1,2) = -1*chi0;
for n=2:NN-1
    H(n,n-1)= -1*chi0;
    H(n,n) = 2*chi0+ V(n);
    H(n,n+1)= -1*chi0;
end
H(NN,NN-1) = -1*chi0;
H(NN,NN) = 2*chi0+V(NN);

% These make it a PBC
%H(1,NN) = -chi0;
%H(NN,1) = -chi0;


[phi,D] = eig(H);
E = diag(D);

% --- Calculate the density matrix from the eigenfunctions
rho = zeros(1,NN);
ff = zeros(NN,1);
for m=1:NN
    ff(m) = 1/(1 + exp((E(m) - mu)/kT)); % Fermi at each eigenenergy.
    for n=1:NN
        rho(n) = rho(n) + ff(m)*phi(n,m)^2;
    end
end

% ------- Calculate the Fermi-Dirac function ----------
Emax = 0.4;
Emin = -0.1;
NE = 250;
EE = zeros(1,NE);
del_E = (Emax-Emin)/NE;

EE = (Emin:del_E:Emax);
fermi = zeros(1,NE);
for m=1:NE
    fermi(m) = 1/(1 + exp((EE(m) - mu)/kT));
end

% --- Calculate the current density via the spectral function
sigma1 = zeros(NN,NN);
sigma2 = zeros(NN,NN);
sig1 = 0.;
sig2 = 0.;
eta = 0;
n = zeros(NN,1);
for m=1:NE
    k = sqrt(2*m0*(EE(m)-V(1))*eV2J)/hbar;
    sig1 = exp(i*k*del_x);
    k = sqrt(2*m0*(EE(m)-V(NN))*eV2J)/hbar;
    sig2 = exp(i*k*del_x);
    sigma1(1,1) = -chi0*sig1;
    sigma2(NN,NN) = -chi0*sig2;
    G = inv( (EE(m) + i*eta)*eye(NN) - H - sigma1 - sigma2);
    n = n + fermi(m)*del_E*real(diag(i*(G-G'))/(2*pi));
end

subplot(2,2,3)
plot(XX,rho,'--k')
hold on
plot(XX,n,'k')
hold off
axis( [ 0 10 0 0.15 ])
set(gca,'fontsize',12)
xlabel('x (nm)')
ylabel('n ')
TT = text(2,.05,'m_1');
set(TT,'FontName','Symbol')
set(TT,'fontsize',12)
TT = text(2.3,.05,sprintf(' = %4.0f meV',1e3*mu));
set(TT,'fontsize',12)
TT = text(2,.02,sprintf('k_BT = %2.0f meV',1e3*kT));
set(TT,'fontsize',12)
grid on

% ffr: this gives warning in Octave
%legend('eigen', 'spec', 2)

