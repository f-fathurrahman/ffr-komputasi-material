clear all
close all

E = linspace(-0.5, 1, 1001);
D = sqrt(E);

hold on
h = plot(D,E,'b');	
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [24])
xlabel('D (E) (arb. units)')
ylabel('E (eV)')
grid on
