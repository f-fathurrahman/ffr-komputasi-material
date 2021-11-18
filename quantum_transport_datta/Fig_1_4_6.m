clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
I0 = q*q/hbar;

%Parameters
U0 = 0.025;
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
dE = E(2)-E(1);

% Lorentzian Density of states per eV
D = (g/(2*pi))./((E.^2) + ((g/2)^2));
D = D./(dE*sum(D)); % Normalizing to one

%Bias
IV = 101;
VV = linspace(0,1,IV);
for iV=1:IV
    Vg = 0;
    Vd = VV(iV);
    %Vd=0; 
    %Vg=VV(iV);
    mu1 = mu;
    mu2 = mu1-Vd;
    UL = -(alphag*Vg)-(alphad*Vd);
    U = 0; %Self-consistent field
    dU = 1;
    while dU > 1e-6
        f1 = 1./(1 + exp((E+ep+UL+U-mu1)./kT));
        f2 = 1./(1 + exp((E+ep+UL+U-mu2)./kT));
        N(iV) = dE*sum(D.*((f1.*g1/g) + (f2.*g2/g)));
        Unew = U0*N(iV);
        dU = abs(U-Unew);
        U = U + 0.1*(Unew-U);
    end
    I(iV) = dE*I0*(sum(D.*(f1-f2)))*(g1*g2/g);
end

hold on
h = plot(VV, N, 'b');
%h = plot(VV, I, 'b' );
xlabel('Voltage (V)')
%ylabel('Current (A)')
ylabel('Number of electrons')
set(h, 'linewidth' ,[2.0])
set(gca, 'Fontsize' ,[25])
grid on