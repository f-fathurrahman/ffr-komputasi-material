% Trans.m.  Calculate the transmision function

clear all

NN = 50;
hbar = 1.06e-34;
m0 = 9.11e-31;
melec = 1.08*m0      % eff. mass of silicon
ecoul = 1.6e-19;
eV2J = 1.6e-19; 
J2eV = 1./eV2J;

del_x = 2.e-10;
DX = del_x*1e9;
X = (DX:DX:NN*DX);
N_center = NN/2;
NC = NN/2;

% Energies are eV
chi0 = J2eV*hbar^2/(2*melec*del_x^2)

%  ---- Specify  the potential ----

V = zeros(1,NN);

% Tunneling barrier
for n=NC-2:NC+2
    V(n) = 0.3;
end

% Resonant barrier
%for n=14:19
for n=NC-7:NC-4
    %V(n) = 0.3;
end 
%for n=31:36
for n=NC+4:NC+7
    %V(n) = 0.3;
end 

subplot(3,2,1)
plot(X,V,'k')
title('Trans')
axis( [ 0 10 -.1 .5 ])
set(gca,'YTick',[ 0 .3 .5 ])
grid on
xlabel(' x (nm)')
ylabel('V (eV)')
set(gca,'fontsize',12)

% ------  Construct the Hamiltonian ---

H = zeros(NN,NN);

H(1,1) = 2*chi0+V(1);
H(1,2) = -1*chi0;

for n=2:NN-1
    H(n,n-1) = -1*chi0;
    H(n,n)   =  2*chi0 + V(n);
    H(n,n+1) = -1*chi0;
end

H(NN,NN-1) = -1*chi0;
H(NN,NN)   =  2*chi0+V(NN);

% --  Specify the energy range ---

Emax = 1;
Emin = 0.;
NE = 250;
EE = zeros(1,NE);
del_E = (Emax-Emin)/NE;
EE = (0:del_E:del_E*(NE-1));

% --- Calculate the transmission function

sigma1 = zeros(NN,NN);
sigma2 = zeros(NN,NN);
gamma1 = zeros(NN,NN);
gamma2 = zeros(NN,NN);
sig1 = 0.;
sig2 = 0.;
eta = 1e-12;
n = zeros(NN,1);
for m=1:NE
    k = sqrt(2*melec*(EE(m)-V(1))*eV2J)/hbar;
    sig1 = exp(i*k*del_x);
    k = sqrt(2*melec*(EE(m)-V(NN))*eV2J)/hbar;
    sig2 = exp(i*k*del_x);
    sigma1(1,1) = -chi0*sig1;
    sigma2(NN,NN) = -chi0*sig2;
    gamma1 = i*(sigma1-sigma1');
    gamma2 = i*(sigma2-sigma2');
    G = inv(  (EE(m) + i*eta)*eye(NN) - H - sigma1 - sigma2);
    TM(m) = real(trace(gamma1*G*gamma2*G'));
end

%subplot(3,2,2)
plot(EE,TM,'k')
grid on
axis( [ 0 1 0 1.2 ])
xlabel('E (eV)')
ylabel('TM')
set(gca,'fontsize',12)
saveas(gcf, 'IMG_trans.png')
