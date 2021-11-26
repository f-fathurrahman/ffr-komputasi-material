clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
eps0 = 8.854E-12;
epsr = 4;
m = 0.25*9.11e-31; %Effective mass
I0 = q*q/hbar;

% Parameters
W = 1e-6;   % width
L = 10e-9;  % length of active region
t = 1.5e-9; % oxide thickness

Cg = epsr*eps0*W*L/t;
Cs = 0.05*Cg;
Cd = 0.05*Cg;
CE = Cg+Cs+Cd;
U0 = q/CE;

alphag = Cg/CE;
alphad = Cd/CE;
%alphag=1; alphad=0.5; U0=0.25;

kT = 0.025;
mu = 0;
ep = 0.2;
v = 1e5;%Escape velocity
g1 = hbar*v/(q*L);
g2 = g1;
g = g1 + g2;
%g1=0.005;g2=0.005;g=g1+g2;

%Energy grid
NE = 501;
E = linspace(-1,1,NE);
dE = E(2) - E(1);

D0 = m*q*W*L/(pi*hbar*hbar);% Step Density of states per eV
D = D0*[zeros(1,251) ones(1,250)];
%D=(2*g/(2*pi))./((E.ˆ2)+((g/2)ˆ2));% Lorentzian Density of states per eV
%D=D./(dE*sum(D));%Normalizing to one

%Reference number of electrons
f0 = 1./(1+exp((E+ep-mu)./kT));
N0 = 2*dE*sum(D.*f0);
ns = N0/(L*W*1e4); %/cmˆ2

% Bias
IV = 61;
VV = linspace(0,0.6,IV);
for iV=1:IV
    Vg = 0.5;
    Vd = VV(iV);
    %Vd=0.5;Vg=VV(iV);
    mu1 = mu;
    mu2 = mu1 - Vd;
    UL = -(alphag*Vg)-(alphad*Vd);
    U = 0; %Self-consistent field
    dU = 1;
    while dU > 1e-6
        f1 = 1./(1+exp((E+UL+U+ep-mu1)./kT));
        f2 = 1./(1+exp((E+UL+U+ep-mu2)./kT));
        N(iV) = dE*sum(D.*((f1.*g1/g)+(f2.*g2/g)));
        Unew = U0*(N(iV)-N0);
        dU = abs(U-Unew);
        U = U + 0.1*(Unew-U);
    end
    I(iV) = dE*I0*(sum(D.*(f1-f2)))*g1*g2/g;
end


hold on
h = plot(VV,I, 'b');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [12])
xlabel('Voltage (V)')
ylabel('Current (A)')
grid on
