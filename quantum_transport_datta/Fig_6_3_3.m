clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
q = 1.602e-19;
a = 5e-10;
L = 10e-9;
k = 0.5*linspace(-1,1,201)/a;
Ek = (hbar^2)*(k.^2)/(2*0.25*m*q);
EE = linspace(0,0.2,201);

%Subband (1,1)
E1 = 2*(hbar^2)*(pi^2)/(2*0.25*m*q*L^2);
M = ((EE-E1) + abs(EE-E1))./(2*abs(EE-E1));

%Subbands (1,2) and (2,1)
E2 = 5*(hbar^2)*(pi^2)/(2*0.25*m*q*L^2);
M = M + (((EE-E2)+abs(EE-E2))./(abs(EE-E2)));

%Subband (2,2)
E3 = 8*(hbar^2)*(pi^2)/(2*0.25*m*q*L^2);
M = M + (((EE-E3)+abs(EE-E3))./(2*abs(EE-E3)));

hold on
h = plot(k, E1+Ek, 'b'); %Part (a)
h = plot(k, E2+Ek, 'b'); %Part (a)
h = plot(k, E3+Ek, 'b'); %Part (a)
%h=plot(M,EE,'b');%Part (b)
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[24])
xlabel('k (/m)'); %Part (a)
%xlabel('M ( E )');%Part (b)
ylabel('E - Ec (eV)');
axis([-1e9 1e9 0 0.3]);%Part (a)
%axis([0 5 0 0.3]);%Part (b)
grid on
