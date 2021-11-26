clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
epsil = 8.854e-12;
q = 1.602e-19;

a0 = 4*pi*epsil*hbar*hbar/(m*q*q);
E0 = q/(8*pi*epsil*a0);
R0 = 0.05*[1:200];
a = (-2*E0)*(1-(exp(-2*R0).*(1+R0)))./R0;
b = (-2*E0)*exp(-R0).*(1+R0);
s = (1 + R0 + ((R0.^2)/3)).*exp(-R0);
Uee = (2*E0)./sqrt(1+(R0.^2));
UNN = (2*E0)./R0;
EB0 = (a+b)./(1+s);
R = a0*R0;

hold on
h = plot(R,EB0, 'b--');
h = plot(R,Uee, 'bx');
h = plot(R,UNN, 'b');
h = plot(R,(2*EB0) + UNN + Uee,' b+');
set(h,'linewidth' ,[2.0])
set(gca,'Fontsize' ,[25])
grid on
xlabel('R (m)' )
ylabel('Energy (eV)' )
axis([0 4e-10 -25 25])
