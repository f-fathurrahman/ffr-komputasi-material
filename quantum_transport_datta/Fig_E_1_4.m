clear all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
I0 = q*q/hbar;

%Parameters
U0 = 0.025;
kT = 0.025;
mu = 0;
ep = 0.2;

N0 = 0;
g1 = 0.005;
g2 = 0.005;
g = g1 + g2;

alphag = 1;
alphad = 0.5;

% Energy grid
NE = 501;
E = linspace(-1,1,NE);
dE = E(2) - E(1);
g1 = 0.005*(E + abs(E))./(E+E+1e-6);% zero for negative E
g2 = 0.005*ones(1,NE);
g = g1 + g2;
%Bias
IV = 101;
VV = linspace(-0.6, 0.6, IV);
for iV = 1:IV
    Vg=0;Vd=VV(iV);
    %Vd=0;Vg=VV(iV);
    mu1 = mu;
    mu2 = mu1 - Vd;
    UL = -(alphag*Vg) - (alphad*Vd);
    U = 0; %Self-consistent field
    dU = 1;
    while dU > 1e-6
        f1 = 1./(1 + exp((E-mu1)./kT));
        f2 = 1./(1 + exp((E-mu2)./kT));
        D = (g./(2*pi))./(((E - ep - UL - U).^2)+((g./2).^2));
        D = D./(dE*sum(D));
        N(iV) = dE*2*sum(D.*((f1.*g1./g) + (f2.*g2./g)));
        Unew = U0*(N(iV) - N0);
        dU = abs(U - Unew);
        U = U + 0.1*(Unew-U);
    end
    I(iV)=dE*2*I0*(sum(D.*(f1-f2).*g1.*g2./g));
end

hold on
%h = plot(VV, N/2, 'b');%Part (a)
h = plot(VV, I, 'b');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
xlabel('Voltage (V)' )
ylabel('Current (A)' )
%ylabel('Number of electrons')
grid on