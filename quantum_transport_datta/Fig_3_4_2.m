clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
I0 = q*q/hbar;

% Parameters
U0 = 0.5;% U0 is 0.25 for part(a), 0.1 for part (b)

kT = 0.025;
mu = 0;
ep = 0.2;
g1 = 0.005;
g2 = 0.005;
g = g1 + g2;
alphag = 1;
alphad = 0.5;

% Bias
IV = 101;
VV = linspace(0,1,IV);
for iV=1:IV
    Vd = 0;
    Vg = VV(iV);
    mu1 = mu;
    mu2= mu1-Vd;
    UL = -(alphag*Vg)-(alphad*Vd);
    f1 = 1/(1 + exp((ep + UL - mu1)/kT));
    f2 = 1/(1 + exp((ep + UL - mu2)/kT));
    f1U = 1/(1 + exp((ep + UL + U0 - mu1)/kT));
    f2U = 1/(1 + exp((ep + UL + U0 - mu2)/kT));
    P1 = ((g1*f1) + (g2*f2))/(1e-6 + (g1*(1-f1)) + (g2*(1-f2)));
    P2 = P1*((g1*f1U) + (g2*f2U))/(1e-6 + (g1*(1 - f1U)) + (g2*(1 - f2U)));
    P0 = 1/(1 + P1 + P1 + P2);
    P1 = P1*P0;
    P2 = P2*P0;
    p0(iV) = P0;
    p1(iV) = P1;
    p2(iV) = P2;
end

hold on
h = plot(VV,p0, 'bo');
h = plot(VV,p1, 'b');
h = plot(VV,p2, 'bx');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
grid on
xlabel('Gate voltage, VG (volts)')
ylabel('Current (Amperes)')
axis([0 1 0 1])