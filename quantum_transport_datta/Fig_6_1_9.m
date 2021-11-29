clear all
close all

t = 3;
kxa = 0;
kyb = linspace(-pi,pi,101);
E1 = (3*t/2)*sqrt(((kxa*2/3).^2)+(((abs(kyb)-(2*pi/3))*2/sqrt(3)).^2));
%a0=b*2/sqrt(3)=a*2/3;
E2 = t*sqrt( 1 + (4*cos(kyb).*cos(kxa)) + (4*cos(kyb).^2) );
k = kyb./pi;

hold on
h = plot(k, E1, 'b');
h = plot(k, -E1, 'b');
h = plot(k, E2, 'bx');
h = plot(k, -E2, 'bx');
axis([-1 1 -15 15])
set(h, 'linewidth', [1.0])
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [24])
xlabel('kyb/pi')
ylabel('Energy (eV)')
grid on