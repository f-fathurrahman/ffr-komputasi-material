clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
q = 1.602e-19;
a = 5e-10;
L = 10e-9;
k = 0.5*linspace(-1,1,201)/a;
Ek = -(hbar^2)*(k.^2)/(2*0.25*m*q);
EE = linspace(0, -0.2, 201);

%Subband (1,1)
E1 = -2*(hbar^2)*(pi^2)/(2*0.25*m*q*L^2);
M = ((E1-EE) + abs(E1-EE))./(2*abs(E1-EE));

%Subbands (1,2) and (2,1)
E2 = -5*(hbar^2)*(pi^2)/(2*0.25*m*q*L^2);
M = M + (((E2-EE)+abs(E2-EE))./(abs(E2-EE)));

%Subband (2,2)
E3 = -8*(hbar^2)*(pi^2)/(2*0.25*m*q*L^2);
M = M + (((E3-EE) + abs(E3-EE))./(2*abs(E3-EE)));

hold on
%h=plot(k,E1+Ek,'b');%Part (a)
%h=plot(k,E2+Ek,'b');%Part (a)
%h=plot(k,E3+Ek,'b');%Part (a)
h = plot(M, EE, 'b');%Part (b)
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [24])
%xlabel('k (/m)');%Part (a)
xlabel('M (E)');%Part (b)
ylabel('E - Ev (eV)');
%axis([-1e9 1e9 -0.3 0]);%Part (a)
axis([0 5 -0.3 0]);%Part (b)
grid on
