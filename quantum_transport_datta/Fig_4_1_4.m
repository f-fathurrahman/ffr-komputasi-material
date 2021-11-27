clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
epsil = 8.854e-12;
q = 1.602e-19;
a0 = 4*pi*epsil*hbar*hbar/(m*q*q);
E0 = q/(8*pi*epsil*a0);

%Basis
L = 0.074e-9/a0;
s = exp(-L)*(1+L+((L^2)/3));
r = linspace(-2e-10,+2e-10,101);
r0 = r/a0;
psi = sqrt(1/(pi*(a0^3)))*(exp(-abs(r0-(L/2)))+exp(-abs(r0+(L/2))));
n = 2*psi.*conj(psi)./(2*(1 + s));
a = -2*E0*(1-((1 + L)*exp(-2*L)))/L;
b = -2*E0*(1+L)*exp(-L);
EB0 = -E0 + ((a+b)/(1 + s));

[a b s EB0]

hold on
h=plot(r,n,'b');
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[25])
grid on
xlabel('R (m)')
ylabel('Electron density (/m^3)')
axis([-2e-10 2e-10 0 2e30])
