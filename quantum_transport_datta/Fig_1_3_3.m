clear all
close all
E = linspace(-0.5, 0.5, 50001);
dE = E(2) - E(1);
gam = 0.05;

% D = (gam/(2*pi))./((E.^2)+((gam/2)^2));

D = (gam/(2*pi))./(((E-0.25).^2) + ((gam/2)^2)); %Use for Fig.1.5.2
D = D + ((gam/(2*pi))./(((E+0.25).^2)+((gam/2)^2))); %Use for Fig.1.5.2

dE*sum(D)

hold on
h = plot(D,E, 'g');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
xlabel('D (E)')
ylabel('E (eV)')
grid on