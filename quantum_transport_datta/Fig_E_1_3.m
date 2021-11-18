clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
I0 = q*q/hbar;

% Parameters
U0 = 0.025;
kT1 = 0.026;
kT2 = 0.025;
ep = 0.2;
g1 = 0.005;
g2 = 0.005;
g = g1 + g2;
alphag = 1;
alphad = 0.5;

% Energy grid
NE = 501;
E = linspace(-1,1,NE);
dE = E(2) - E(1);

g1 = 0.005*(E + abs(E))./(E + E + 1e-6); % zero for negative E
g2 = 0.005*ones(1,NE);
g1 = g2;
g = g1 + g2;

% Bias
IV = 101;
VV = linspace(-0.25,0.25,IV);
for iV = 1:IV
    mu1 = ep + VV(iV);
    mu2 = mu1;
    f1 = 1./(1+exp((E - mu1)./kT1));
    f2 = 1./(1+exp((E - mu2)./kT2));
    D = (g./(2*pi))./(((E-ep).^2)+((g./2).^2));
    D = D./(dE*sum(D));
    I(iV) = dE*2*I0*(sum(D.*(f1-f2).*g1.*g2./g));
end

hold on
%h=plot(VV,N/2, ' b ' );%Part (a)
h = plot(VV, I, 'b');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
xlabel('Voltage (V)')
ylabel('Current (A)')
%ylabel( ' Number of electrons ---> ' )
grid on
