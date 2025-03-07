% Current.m.  Calculate current for a given voltage

clear all

NN = 50;
NC = NN/2;
hbar = 1.054e-34;
m0 = 9.1e-31;
melec = 1.08*m0;       % Eff. mass of silicon
ecoul = 1.6e-19;
eV2J = 1.6e-19; 
J2eV = 1./eV2J;

del_x = 2.e-10;         % The cell size
DX = del_x*1e9;
XX = (DX:DX:NN*DX);

G0= ecoul^2/(2*pi*hbar);    % Quantum conductance

% Energies are eV
chi0 = J2eV*hbar^2/(2*melec*del_x^2)

mu = input('mu (left side) -->')
kT = 0.025;

%  --- Channel potential --------------

V = zeros(1,NN);

% Tunneling barrier
for n=NC-2:NC+2
    %  V(n) = 0.3;
end

%Resonant barrier
for n=NC-7:NC-4
    %    V(n) = 0.3;
end 
for n=NC+4:NC+7
    %    V(n) = 0.3;
end 

% Triangle potential
for n=15:25
    %V(n) = 0.05*(n-15);
end
for n=25:35
    %V(n) = 0.5*(1-.1*(n-25));
end

% --------------------------------------

VDS = input('Voltage across the channel -->')

for n=1:NN
    VD(n) = -(n-1)*VDS/(NN-1);
    %VD(n) = -.5*(n-NC)*VDS/(NN-1);
    V(n) = V(n) + VD(n);
end

Vmin = min(V);

subplot(3,2,1)
plot(XX,V,'k')
title('Current')
axis( [ 0 10 Vmin .6 ])
grid on
TT = text(3,.4,sprintf('V_D_S = %5.2f V',VDS));
set(TT,'fontsize',12)
set(gca,'fontsize',12)
xlabel(' x (nm)')
ylabel('V (eV)')

% -- Construct the Hamiltonian --

H = zeros(NN,NN);

H(1,1) = 2*chi0+V(1);
H(1,2) = -1*chi0;

for n=2:NN-1
    H(n,n-1)= -1*chi0;
    H(n,n)  =  2*chi0+ V(n);
    H(n,n+1)= -1*chi0;
end

H(NN,NN-1) = -1*chi0;
H(NN,NN)   =  2*chi0+V(NN);

% -- Calculate the Fermi functions at the contacts ---

Emax = 1;
Emin = 0.;
NE = 250;
EE = zeros(1,NE);
del_E = (Emax-Emin)/NE;
EE = (0:del_E:del_E*(NE-1));

fermi1 =zeros(1,NE);
fermi2 =zeros(1,NE);
TM =zeros(1,NE);
over =zeros(1,NE);
for m=1:NE
    fermi1(m) = 1/(1 + exp( (EE(m) - (mu + V(1))) /kT ));
    fermi2(m) = 1/(1 + exp( (EE(m) - (mu + V(NN))) /kT ));
end

subplot(3,2,2)
plot(EE,fermi1,'k')
hold on
plot(EE,fermi2,'k--')
hold off
set(gca,'fontsize',12)
ylabel('Fermi')
xlabel('E (eV)')
legend('f_1','f_2')
TT = text(.1,.6,'m_1','FontName','Symbol');
set(TT,'fontsize',12)
TT = text(.15,.61,sprintf(' = %4.3f eV',mu+V(1)));
set(TT,'fontsize',12)
TT = text(.1,.3,'m_2','FontName','Symbol');
set(TT,'fontsize',12)
TT = text(.15,.36,sprintf(' = %4.3f eV',mu+V(NN)));
set(TT,'fontsize',12)
%TT = text(.1,.1,sprintf('kT = %4.3f eV',kT));
%set(TT,'fontsize',12)
grid on

% --- Calculate the Transmission function and the current

sigma1 = zeros(NN,NN);
sigma2 = zeros(NN,NN);
gamma1 = zeros(NN,NN);
gamma2 = zeros(NN,NN);
sig1 = 0.;
sig2 = 0.;
eta = 0;
n = zeros(NN,1);
I = 0.;
for m=1:NE
    k = sqrt(2*melec*(EE(m)-V(1))*eV2J)/hbar;
    sig1 = exp(i*k*del_x);
    k = sqrt(2*melec*(EE(m)-V(NN))*eV2J)/hbar;
    sig2 = exp(i*k*del_x);
    sigma1(1,1) = -chi0*sig1;
    sigma2(NN,NN) = -chi0*sig2;
    gamma1 = i*(sigma1-sigma1');
    gamma2 = i*(sigma2-sigma2');
    G = inv(  (EE(m) + i*eta)*eye(NN) - H - sigma1 - sigma2);
    TM(m) = real(trace(gamma1*G*gamma2*G'));
    I = I + (eV2J*del_E)*(ecoul/(2*pi*hbar))*TM(m)*(fermi1(m) - fermi2(m));
    %I = I + G0*del_E*TM(m)*(fermi1(m) - fermi2(m));
    over(m) = TM(m)*(fermi1(m) - fermi2(m));
end
I

subplot(3,2,3)
plot(EE,TM,'k')
grid on
axis( [ 0 1 0 1.2 ])
set(gca,'fontsize',12)
ylabel('TM')
xlabel('E (eV)')
%TT = text(.4,.3,sprintf('V_D_S = %5.2f V',VDS))'
%set(TT,'fontsize',12)
%TT = text(.4,.1,sprintf('I = %5.3f uA',1e6*I))'
%set(TT,'fontsize',12)

subplot(3,2,4)
plot(EE,fermi1-fermi2,'k--')
hold on
%plot(EE,over,'k')
bar(EE,over)
hold off
%legend('f_1 - f_2','T*(f_1-f_2)')
legend('f_1 - f_2')
%TT = text(.5,.75,sprintf('F_t_o_t = %5.3f ',G0*sum(fermi1-fermi2)));
%set(TT,'fontsize',12)
TT = text(.5,.4,sprintf(' I  = %5.2f uA',1e6*I));
set(TT,'fontsize',12)
TT = text(.5,.2,sprintf(' G = %5.2f uS',1e6*(I/VDS)));
set(TT,'fontsize',12)
ylabel('(f_1 - f_2) x TM')
set(gca,'fontsize',12)
xlabel('E (eV)')
axis([ 0 1 0 1 ])
grid on
%saveas(gcf,'trans.png')
