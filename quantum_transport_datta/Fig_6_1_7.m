clear all
close all

t = 3;
m = 65; %66 for (a), 65 for (b)
D = 2*m*0.14*sqrt(3)/(2*pi);
Eg = 2*t*0.14/D;
nu = round(2*m/3) + 0; % +1 is used for higher mode
kyb = 2*pi*nu/(2*m);
kxa = 0.05*linspace(-pi,pi,101);
E1 = (3*t/2)*sqrt(((kxa*2/3).^2)+(((abs(kyb)-(2*pi/3))*2/sqrt(3)).^2));
%a0=b*2/sqrt(3)=a*2/3;
E2 = t*sqrt(1+(4*cos(kyb).*cos(kxa))+(4*cos(kyb).^2));
k = kxa./pi;

[D Eg nu min(E1)]

hold on
h = plot(k, E1, 'b');
h = plot(k, -E1, 'b');
axis([-0.05 0.05 -0.6 0.6])
set(h, 'linewidth', [1.0])
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [24])
xlabel('kxa/pi (fraction of maximum value')
ylabel('Energy (eV)')
grid on
