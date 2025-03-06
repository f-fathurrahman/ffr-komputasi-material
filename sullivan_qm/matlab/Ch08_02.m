% Spectral.m This program calculates the spectral function
clear all
hbar = 1.054e-34;
m0 = 9.11e-31;
ecoul = 1.6e-19;
eV2J = 1.6e-19;
J2eV = 1./eV2J;
del_x = 1.e-10;
NN = 100;
V = zeros(1,NN);

% Energies are in eV
chi0 = J2eV*hbar^2/(2*m0*del_x^2);

% Convert to nm for plotting
DX = del_x*1e9;
XX = (DX:DX:NN*DX);

% ------- Construct the Hamiltonian ----------
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

% ---------------------------------------------
% This MATLAB command calculates the eigenenergies
% and correspondig eigenfunctions
[phi,D] = eig(H);
E = diag(D);

% ------- Calculate the broadened DOS ---------
Ein = 0.03375;
gamma = 0.0001;

wtot = 0.;
del_E = 0.0005;
EE = (del_E:del_E:del_E*100);
func = zeros(1,100);
for m=1:6
    for l=1:100
        func(l) = func(l) + del_E*(gamma/(2*pi))/( ( gamma/2)^2 + (E(m) - EE(l))^2);
        wtot = wtot + func(l);
    end
end
wtot

% -------- Calculate the spectral function -----------------
AA = zeros(NN,NN);
for m=1:5
    wide = del_E*(gamma/(2*pi))/( ( gamma/2)^2 + (Ein - E(m))^2);
    for n=1:NN
        for l=1:NN
            AA(n,l) = AA(n,l) + wide*phi(n,m)*phi(l,m);
        end
    end
end

subplot(2,2,1)
mesh(XX,XX,AA)
ylabel('x^1 (nm)')
xlabel('x (nm)')
zlabel('A(E)')
view(-20,30)
TT = text(5,8,.0006,sprintf('E = %6.4f eV',Ein))
set(TT,'fontsize',12)
set(gca,'fontsize',12)
title('Spectral')