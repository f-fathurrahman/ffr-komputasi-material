% Spec.m. Calculate the electron density in the 10 nm inﬁnite well
% ﬁrst by the eigenfunction, and then with the
% spectral matrix.

hbar = 1.06e-34;
m0 = 9.1e-31;
elec = 1.6e-19;
eV2J = 1.6e-19;
J2eV = 1./eV2J;

mu = 0.01; % ffr: also try 0.04 and 0.1
kT = 0.001;
%kT = 0.0259; % Room temp

% del_x is the cell size (m)
del_x = 1.e-10;
NN = 100;
N_center = NN/2;

% Energies are eV
chi0 = (1/elec)*hbar^2/(2*m0*del_x^2);

% Convert to nm for plotting
DX = del_x*1e9;
XX = (DX:DX:NN*DX);
mid = DX*N_center;

% -------- Specify the potential---------
V = zeros(1,NN);
for n=1:NN
    %V(n) = .01*((N_center-n)/NN);
end

% ------ Construct the Hamiltonian ------
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

% Solve eigenvalue equation
[phi,D] = eig(H);
E = diag(D);

% --- Calculate the density matrix from the eigenfunctions
rho = zeros(1,NN);
ff = zeros(NN,1);
for m=1:100
    ff(m) = 1/(1 + exp((E(m) - mu)/kT)); % Fermi at each eigenenergy.
    for n=1:NN
        rho(n) = rho(n) + ff(m)*phi(n,m)^2;
    end
end

% -- Calculate the Fermi-Dirac function ----------
Emax = 0.4;
Emin = 0.00001;
NE = 250;
EE = zeros(1,NE);
del_E = (Emax-Emin)/NE;
EE = (del_E:del_E:del_E*NE);
fermi =zeros(1,NE);
for m = 1:NE
    fermi(m) = 1/(1 + exp((EE(m) - mu)/kT));
end

% --- Calculate the electron density via the spectral function
sigma = zeros(NN,NN);
eta = 2e-3;
% This is an artiﬁcial loss term
n = zeros(NN,1);
for m = 1:100
    k = sqrt(2*m0*EE(m)*eV2J)/hbar;
    sig1 = exp(i*k*del_x);
    %sig1 = 0.;
    % sig1 = 0. means no contact potential
    sigma(1,1) = -chi0*sig1;
    sigma(NN,NN) = -chi0*sig1;
    G = inv( (EE(m) + i*eta)*eye(NN) - H - sigma);
    n = n + fermi(m)*del_E*real(diag(i*(G-G'))/(2*pi));
end


amax = 0.1;
%subplot(3,2,1)
plot(XX,rho,'--k')
hold on
plot(XX,n,'k')
hold off
axis( [ 0 10 0 amax ])
set(gca,'fontsize',12)
xlabel('x (nm)')
ylabel('n ')
TT = text(1,.7*amax,'m','FontName','Symbol');
set(TT,'fontsize',12)
TT = text(1.2,.7*amax,sprintf(' = %3.0f meV',1e3*mu));
set(TT,'fontsize',12)
legend('eigen','spec')
title('Spec')

