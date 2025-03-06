% Cal_rho.m Calculate the electron density and the current density.
clear all
hbar = 1.054e-34;
m0 = 9.11e-31;
ecoul = 1.6e-19;
eV2J = 1.6e-19;
J2eV = 1./eV2J;
hbar_eV = J2eV*hbar
del_x = 1.e-10;
NN = 100;
N_center = NN/2;

% Energies are meV
chi0 = J2eV*hbar^2/(2*m0*del_x^2);


% Convert to nm for plotting
DX = del_x*1e9;
X = (DX:DX:NN*DX);
% This adds the H.O. Potential
U = zeros(1,NN);
Eref = .1;
% Ground state energy in eV
k = m0*(Eref/hbar_eV)^2;
for n=1:NN
    % U(n) = J2eV*.5*k*((N_center-n)*a)^2;
    % ffr: if commented this means that there is no external pot?
end

% ---- Construct the Hamiltonian ---
H = zeros(NN,NN);
H(1,1) = 2*chi0 + U(1);
H(1,2) = -1*chi0;
for n=2:NN-1
    H(n,n-1) = -1*chi0;
    H(n,n) = 2*chi0+ U(n);
    H(n,n+1) = -1*chi0;
end
H(NN,NN-1) = -1*chi0;
H(NN,NN) = 2*chi0+U(NN);
% ----------------------------------
% This MATLAB command calculated the eigenenergies
% and corresponding eigenfunctions
[phi,D] = eig(H);
E = diag(D);

% ------- Calculate the Fermi function ------
f = zeros(1,NN);
mu = 0.05;
% Chemical potential
kT = 0.001;
for m=1:NN
    f(m) = 1./(1+exp((E(m)-mu)/kT));
end

subplot(4,2,1)
plot(E,f,'o--k')
axis( [ 0 .1 0 1.2])
xlabel(' E (eV)')
ylabel('f(E)')
grid on
title('Cal-rho')
TT = text( .05, .8, 'm ','FontName','Symbol' );
set(TT,'fontsize',12)
TT = text( .053, .8, [' = ',num2str(mu) ]);
set(TT,'fontsize',12)
TT = text( .05, .3, ['k_BT = ',num2str(kT) ]);
set(TT,'fontsize',12)
set(gca,'fontsize',12)

% ------- Calculate rho ----------
rho = zeros(1,NN);
rho2 = zeros(NN,NN);
for m=1:20
    for n=1:NN
        rho(n) = rho(n) + f(m)*phi(n,m)^2;
        for l=1:NN
            rho2(l,n) = rho2(l,n) + f(m)*phi(n,m)*phi(l,m);
        end
    end
end

subplot(4,2,3)
plot(X,rho,'k')
xlabel(' x (nm)')
ylabel('n')
set(gca,'fontsize',12)
subplot(2,2,3)
mesh(X,X,rho2)
view(-15,30)
axis( [ 0 10 0 10 -0.005 0.04 ])
zlabel('r','FontName','Symbol')
xlabel('x (nm)')
ylabel('x^l(nm)')
set(gca,'fontsize',12)


