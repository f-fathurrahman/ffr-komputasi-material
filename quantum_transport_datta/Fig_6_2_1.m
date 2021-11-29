clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
q = 1.602e-19;
L = 1e-9;

D2 = zeros(1,101);
Lz=20e-9;%5e-9 for (a),20e-9 for (b)
E0=(hbar^2)*(pi^2)/(2*q*m*Lz^2);
for p = 1:25
    E = linspace(-0.1,0.25,101);
    thet = (E + abs(E))./(2*E);
    EE = E-(p*p*E0);
    theta = (EE+abs(EE))./(2*EE);
    D1 = (L)*q*m*thet.*real((2*m*E*q).^(-0.5))./(pi*hbar);
    D2 = D2 + ((L^2)*q*m*theta./(2*pi*hbar*hbar));
    D3 = (L^3)*q*m*thet.*real((2*m*E*q).^0.5)./(2*pi*pi*hbar*hbar*hbar);
end

hold on
h = plot(D2,E,'b');
h = plot(D3.*Lz/L,E,'b');
%axis([0 10 -0.1 0.25]);%Part (a)
axis([0 40 -0.1 0.25]);%Part (b)
set(h,'linewidth',[1.0])
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[24])
xlabel('D(E) (per eV per nm^2)')
ylabel('Energy (eV)')
grid on
