%E.3.5c: Unrestricted scf
clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
I0 = q*q/hbar;

% Parameters
U0 = 0.25;
kT = 0.025;
mu = 0;
ep = 0.2;
g1 = 0.005;
g2 = 0.005;
g = g1 + g2;
alphag = 1;
alphad = 0.5;

%Energy grid
NE = 501;
E = linspace(-1,1,NE);
dE = E(2) - E(1);
D = (g/(2*pi))./((E.^2) + ((g/2)^2)); % Lorentzian Density of states per eV
D = D./(dE*sum(D)); % Normalizing to one

% Bias
IV = 101;
VV = linspace(0,1,IV);
for iV=1:IV
    Vg = 0;
    Vd = VV(iV);
    %Vd=0;Vg=VV(iV);
    mu1 = mu;
    mu2 = mu1 - Vd;
    UL = -(alphag*Vg)-(alphad*Vd);
    Uup = 0;
    Udn = 0.1;%Unrestricted self-consistent field
    dU = 1;
    while dU > 0.001
        f1up = 1./(1 + exp((E + ep + UL + Uup - mu1)./kT));
        f2up = 1./(1 + exp((E + ep + UL + Uup - mu2)./kT));
        Nup(iV) = dE*sum(D.*((f1up.*g1) + (f2up.*g2))./(g1+g2));
        f1dn = 1./(1 + exp((E + ep + UL + Udn - mu1)./kT));
        f2dn = 1./(1 + exp((E + ep + UL + Udn - mu2)./kT));
        Ndn(iV) = dE*sum(D.*((f1dn.*g1) + (f2dn.*g2))./(g1+g2));
        Udnnew = 2*U0*(Nup(iV) - 0.5);
        Udn = Udn + 0.1*(Udnnew - Udn);
        Uupnew = 2*U0*(Ndn(iV) - 0.5);
        Uup = Uup + 0.1*(Uupnew - Uup);
        dU = abs(Uup - Uupnew) + abs(Udn - Udnnew);
    end
    Iup(iV) = dE*I0*sum(D.*(f1up - f2up))*(g1*g2/g);
    Idn(iV) = dE*I0*sum(D.*(f1dn - f2dn))*(g1*g2/g);
end

hold on

%h = plot(VV, Iup + Idn, 'b');
h=plot(VV,Nup,'bo');%Part (b)
h=plot(VV,Ndn,'bx');%Part (b)

set(h, 'linewidth', [2.0])
set(gca,'Fontsize', [25])
xlabel('Voltage (V)')
%ylabel('Current (A)')
ylabel('Number of electrons'); % Part (b)

grid on
