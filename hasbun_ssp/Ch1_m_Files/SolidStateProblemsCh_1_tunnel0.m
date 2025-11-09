%copyright by J. E Hasbun and T. Datta
%tunnel0.m
%Here we use the approximate tunneling formula with maximum
%amplitude of unity and plot it versus energy and spacing.
clear
hbC=197.0;    %hbar*C in eV nm
V0=4.3;       %example work function in eV
% ****** varying the applied voltage ******************
Va=0:V0/100.0:V0;  %applied voltage in electron volts
E=Va;                            %energy in eV
mc2=0.511e6;                     %rest mass of the electon mc^2 in eV
k=sqrt(2*mc2*(V0-E)/hbC^2);
d=0.3;                           %gap width in nm
T1_approx=exp(-2.0*k*d);         %approximation
subplot(1,2,1)
plot(E,T1_approx,'k--','LineWidth',2)
axis([0 1.1*V0 0 1])
text(0.1,0.8,'T(E)=e^{(-2.0*k(E)*d)}','FontSize',16)
xlabel('E (eV)','FontSize',16), ylabel ('T','FontSize',16)
title('T vs E','FontSize',16)
% ****** varying the gap width ******************
s=0:d/25:d;
EE=0.5*V0;
k=sqrt(2*mc2*(V0-EE)/hbC^2);
T2_approx=exp(-2*k*s);           %approximation
subplot(1,2,2)
plot(s,T2_approx,'k--','LineWidth',2)
axis([0 1.1*d 0 1])
text(d/4,0.4,'T(d)=e^{(-2.0*k*d)}','FontSize',16)
xlabel('d (nm)','FontSize',16), ylabel ('T','FontSize',16)
title('T vs a','FontSize',16)
