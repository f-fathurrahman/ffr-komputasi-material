clear all
close all

k = linspace(-1,1,21);
a = 2;
b = 1;
E1 = sqrt( a^2 + b^2 + 2*a*b.*cos(pi*k) );
hold on
h = plot(k, E1, 'b');
h = plot(k, -E1, 'b');
set(h,'linewidth', [2.0])
set(gca,'Fontsize', [25])
xlabel('k (in units of pi/a)')
ylabel('Energy (eV)')
grid on