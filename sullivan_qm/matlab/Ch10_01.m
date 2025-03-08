%find_var.m  Estimate an eigenfunction
%            through the variational method.

clear all

NN = 100;
hbar = 1.054e-34;
melec = 9.11e-31;
ecoul = 1.6e-19;
eV2J = 1.6e-19; 
J2eV = 1./eV2J;

del_x = .1e-9;
DX = del_x*1e9;
XX = (DX:DX:NN*DX);
mid = DX*NN/2;

% Energy is in J.
chi0 = hbar^2/(2*melec*del_x^2);

%  ---- Potential -----
V = zeros(1,NN);

% V-shapted potential
Eref  = .05*eV2J; 
for n=1:NN
    V(n) = (Eref/(NN/2))*abs(n-NN/2);
end

% Harmonic oscillator
Eref = 0.02*eV2J;
for n=1:NN
    %V(n) = 0.5*melec*(Eref/hbar)^2*del_x^2*(n-NN/2)^2;
end

subplot(3,2,1)
plot(XX,J2eV*V,'k');
set(gca,'fontsize',12)
ylabel('V (eV)')
xlabel('nm')
title('Find-var')
Vmax = max(J2eV*V);
grid on
axis( [ 0  DX*NN 0 .06 ])

% ------- Create the Hamiltonian matrix ---

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
% These two lines add the PBC
%H(1,NN) = -1*chi0
%H(NN,1) = -1*chi0

% ------------------------------------
     
% This calculates the eigenfunction and
%   eigenenergies
[phi,D] = eig(H);    

% Write the eigevalues in meV.
for m = 1:NN
E(m) = 1e3*J2eV*D(m,m);
end

% -----  Plot an eigenfunction ----------
subplot(3,2,2)
plot(XX,phi(:,2),'k')
TT = text( 3.5,.075,sprintf('%5.2f meV',E(2)));
set(TT,'fontsize',12);
TT = ylabel('f_2','FontName','Symbol','fontsize',12);
axis( [ 0 DX*NN -.25 .25 ])
xlabel('nm')
grid on
set(gca,'fontsize',12)

%  --- Guess at the eigenfunction  -------

prl = zeros(1,NN);

nc = NN/2;
sigma = 1.;
while sigma > 0
    sigma = input('Sigma -->')
    LL = input('LL(nm) -->')
    LL = LL/DX;
    for n=2:NN-1
        %prl(n) = exp(-.5*((n-nc)/sigma)^2);       % Ground state
        prl(n) = exp(-.5*((n-nc)/sigma)^2)*sin(2*pi*(n-nc)/LL);  % 2nd state
    end
    ptot = prl*prl';
    prl = prl/sqrt(ptot);
    prl*prl';

    subplot(3,2,5)
    plot(XX,prl,'k')
    TT = text(4,.2,'s','FontName','Symbol')
    set(TT,'fontsize',12);
    TT = text( 4.5,.2,sprintf(' =  %4.2f ',sigma));
    set(TT,'fontsize',12);
    %TT = text( 4.,.0,sprintf('L  =  %4.2f ',LL*0.1));
    %set(TT,'fontsize',12);

    Evar = J2eV*prl*H*prl'
    TT = text( 2.5,-.2,sprintf('E_v_a_r = %5.2f meV',1e3*Evar));
    set(TT,'fontsize',12);
    axis( [ 0 DX*NN -0.25 .25 ])
    xlabel('nm')
    grid on
    set(gca,'fontsize',12)
    title('Find-var')

    saveas(gcf,'var.png')
end
