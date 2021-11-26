clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
I0 = q*q/hbar;

% Parameters
U0 = 0.1; % U0 is 0.25 for part(a), 0.025 for part (b)
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
VV = linspace(0, 1.5, IV);
for iV = 1:IV
    Vg = 0;
    Vd = VV(iV);
    %Vd = 0; Vg=VV(iV);
    mu1 = mu;
    mu2 = mu1 - Vd;
    UL = -(alphag*Vg)-(alphad*Vd);
    % Multielectron method
    f1  = 1/(1 + exp((ep + UL - U0/2 - mu1)/kT));
    f2  = 1/(1 + exp((ep + UL - U0/2 - mu2)/kT));
    f1U = 1/(1 + exp((ep + UL + U0/2 - mu1)/kT));
    f2U = 1/(1 + exp((ep + UL + U0/2 - mu2)/kT));
    P1 = ((g1*f1) + (g2*f2))/(1e-6 + (g1*(1 - f1)) + (g2*(1 - f2)));
    P2 = P1*((g1*f1U)+(g2*f2U))/(1e-6 + (g1*(1 - f1U)) + (g2*(1 - f2U)));
    P0 = 1/(1 + P1 + P1 + P2);
    P1 = P1*P0;
    P2 = P2*P0;
    I1(iV) = 2*I0*((P0*g1*f1) - (P1*g1*(1 - f1)) + (P1*g1*f1U) - (P2*g1*(1 - f1U)));
    I2(iV) = 2*I0*((P0*g2*f2) - (P1*g2*(1 - f2)) + (P1*g2*f2U) - (P2*g2*(1 - f2U)));
end

%RSCF method (same as Fig.1.4.6 with added factor of two)

%Energy grid
NE = 501;
E = linspace(-1,1,NE);dE=E(2)-E(1);
D = (g/(2*pi))./((E.^2)+((g/2)^2)); % Lorentzian Density of states per eV
D = D./(dE*sum(D)); %Normalizing to one

%Bias
for iV=1:IV
    Vg=0;Vd=VV(iV);
    %Vd=0;Vg=VV(iV);
    mu1 = mu;
    mu2 = mu1 - Vd;
    UL = -(alphag*Vg) - (alphad*Vd);
    U = 0; %Self-consistent field
    dU = 1;
    while dU > 1e-6
        F1 = 1./(1 + exp((E + ep + UL + U - mu1)./kT));
        F2 = 1./(1 + exp((E + ep + UL + U - mu2)./kT));
        N(iV) = dE*2*sum(D.*((F1.*g1/g) + (F2.*g2/g)));
        Unew = U0*N(iV);
        dU = abs(U - Unew);
        U = U + 0.1*(Unew-U);
    end
    I(iV) = dE*2*I0*(sum(D.*(F1-F2)))*(g1*g2/g);
end

hold on
h = plot(VV, I1, 'b');
h = plot(VV, I, 'b--');
set(h,'linewidth', [2.0])
set(gca,'Fontsize', [25])
grid on
xlabel('Drain Voltage, VD (volts)')
ylabel('Current (Amperes)')
axis([0 1.5 0 1.4e-6])
