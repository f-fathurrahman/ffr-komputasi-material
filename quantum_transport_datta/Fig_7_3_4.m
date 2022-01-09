clear all
close all

E = linspace(-.25,.25,501);
dE = E(2) - E(1);
kT = 0.025;
Ef = 0;
V = 0;
mu1 = Ef + (V/2);
mu2 = Ef - (V/2);
f1 = 1./(1+exp((E-mu1)./kT));
f2 = 1./(1+exp((E-mu2)./kT));
FT = [0 diff(f1)];
FT = FT.*(-1/dE);
%dE*(sum(f1-f2))/V

hold on
h = plot(FT, E, 'b');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [24])
grid on
