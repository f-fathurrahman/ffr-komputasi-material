% Broadening.m  Calculate the thermal broadening function

clear all

hbar = 1.06e-34;
m0 = 9.1e-31;
melec = 1.08*m0;
elec = 1.6e-19;
eV2J = 1.6e-19; 
J2eV = 1./eV2J;

mu = 0.25;
kT = 0.01;

% --- Calculate the Fermi-Dirac and the broadening functions ------

Emax = .4;
Emin = 0;
NE = 250;
EE = zeros(1,NE);
del_E = (Emax-Emin)/NE;
EE = (0:del_E:del_E*(NE-1));

fermi =zeros(1,NE);
derf =zeros(1,NE);
fermi(1) = 1;
for m=2:NE
    fermi(m) = 1/(1 + exp((EE(m) - mu)/kT));
    derf(m) = -(fermi(m) - fermi(m-1) )/del_E;
end

subplot(3,2,1)
plot(EE,fermi,'k')
title('Broadening')

TT = text(.05,.7,'m ','FontName','Symbol');
set(TT,'fontsize',12)
TT = text(.07,.7,sprintf(' = %5.3f eV',mu))
set(TT,'fontsize',12)
TT = text(.05,.3,sprintf('kT = %5.3f eV',kT))
set(TT,'fontsize',12)
ylabel('f_T (E)')
set(gca,'fontsize',12)
grid on

dmax = max(derf);
subplot(3,2,3)
plot(EE,derf,'k')
axis( [ 0 .4 0 dmax])
%saveas(gcf,'rho.bmp')
ylabel('F_T(E)')
xlabel('E (eV)')
ylabel('F_T (E)')
set(gca,'fontsize',12)
grid on

%saveas(gcf,'rho.png')
