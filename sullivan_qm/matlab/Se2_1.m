clear all; close all;

NN = 800;
hbar = 1.054e-34;
melec = 9.1e-31;
eV2J = 1.6e-19;
J2eV = 1/eV2J;

del_x = 0.05e-9; % The cells size
dt = 0.5e-17; % Time steps
Ra = (0.5*hbar/melec)*(dt/del_x^2) % ra must be < .1
DX = del_x *1e9; % Cell size in nm.
XX = (DX:DX:DX *NN); % Length in nm for plotting

chi0 = hbar^2/(2*melec*del_x^2);

V = zeros(1,NN);

% V shaped potential
for n=1:NN
    %V(n) = eV2J *(.0005)*(abs(NN/2-n));
end

% E ﬁeld
for n = 1:NN
    %V(n) = eV2J *(.002)*(NN/2-n);
end

subplot(3,2,1)
plot(XX,J2eV*V,'k');
set(gca,'fontsize',12)
ylabel('V (eV)')
xlabel('nm')
Umax = max(J2eV *V);
title('Se2-1')
%axis( [ 0 DX *NN 0 Umax ])

% ----- Eigenvalue calculation -------------------------------
% Specify the Hamiltonian
H = zeros(NN,NN);
H(1,1) = 2 *chi0 + V(1);
H(1,2) = -1*chi0;
for n = 2:NN-1
    H(n,n-1) = -1*chi0;
    H(n,n) = 2*chi0 + V(n);
    H(n,n+1)= -1*chi0;
end
H(NN,NN-1) = -1*chi0;
H(NN,NN) = 2*chi0+V(NN);

% Switch to PBC
%H(NN,1) = -1*chi0;
%H(1,NN) = -1*chi0;

[phi,D] = eig(H);

%Plot the eigenfunctions
subplot(3,2,3)
plot(XX,phi(:,1),'k')
TT = ylabel( 'f_1','FontName','Symbol','fontsize',12)
TT = text(5,0.05,sprintf( '%7.4f eV ',J2eV*D(1,1)));
set(TT,'fontsize',12)
set(gca,'fontsize',12)
title('Se2-1')
subplot(3,2,4)
plot(XX,phi(:,2),'k')
TT = ylabel('f_2','FontName','Symbol','fontsize',12)
TT = text(5, 0.03,sprintf('%7.4f eV ',J2eV*D(2,2)));
set(TT,'fontsize',12)
set(gca,'fontsize',12)
subplot(3,2,5)
plot(XX,phi(:,3),'k')
TT = ylabel('f_3', 'FontName', 'Symbol', 'fontsize',12)
TT = text(5, 0.03, sprintf('%7.4f eV',J2eV*D(3,3)));
set(TT,'fontsize',12)
set(gca,'fontsize',12)
xlabel('nm')
subplot(3,2,6)
plot(XX,phi(:,4),'k')
TT = ylabel( 'f_4','FontName','Symbol','fontsize',12)
TT = text(5,.03,sprintf( '%7.4f eV ',J2eV*D(4,4)));
set(TT,'fontsize',12)
set(gca,'fontsize',12)
xlabel('nm')

% Initialize the FDTD simulation
sigma = 40; % Pulse width
lambda = 40; % Pulse wavelength
nc = NN/2 - 30; % Starting position of the pulse
prl = zeros(1,NN);
pim = zeros(1,NN);
ptot = 0.0;
for n = 2:NN-1
    %prl(n) = exp( -1.*((n-nc)/sigma)^2)*cos(2*pi*(n-nc)/lambda) ;
    %pim(n) = exp( -1.*((n-nc)/sigma)^2)*sin(2*pi*(n-nc)/lambda) ;
    prl(n) = phi(n,1); % This initializes an eigenstate
    ptot = ptot + prl(n)^2 + pim(n) ^2;
end
pnorm = sqrt(ptot);

% Normalize and check
ptot = 0.0;
for n=1:NN
    prl(n) = prl(n)/pnorm;
    pim(n) = pim(n)/pnorm;
    ptot = ptot + prl(n)^2 + pim(n)^2;
end

ptot

T = 0;
n_step = 2000;

for m = 1:n_step
    T
    T = T + 1;
    for n = 2:NN-1
        prl(n) = prl(n) - Ra *(pim(n-1) -2*pim(n) + pim(n +1)) + (dt/hbar) *V(n)*pim(n);
    end
    for n = 2:NN-1
        pim(n) = pim(n) + Ra *(prl(n-1) -2*prl(n) + prl(n +1)) - (dt/hbar) *V(n)*prl(n);
    end
end


% Check normalization
ptot1 = 0.;
for n =1:NN
    ptot1 = ptot1 + prl(n)^2 + pim(n)^2;
end
ptot1

% Calculate the expected values
PE = 0.;
for n = 1:NN
    psi(n) = prl(n) + i *pim(n);
    PE = PE + psi(n)*psi(n)' * V(n);
end
PE = PE*J2eV;

ke = 0.0 + j*0.0;
for n = 2:NN-1
    lap_p = psi(n +1) - 2 *psi(n) + psi(n -1);
    ke = ke + lap_p *psi(n)';
end
KE = -J2eV*((hbar/del_x)^2/(2*melec))*real(ke);

subplot(1,1,1) % This creates a new window
subplot(3,2,1)
plot(XX,prl,'k')
hold on
plot(XX,pim,'--k')
plot(XX,J2eV*V,'-.k')
hold off
axis( [ 0 DX *NN -.25 .30 ])
TT = text(1,.2,sprintf( '%5.3f ps',T*dt*1e12));
set(TT,'fontsize',12)
TT = text(1, -.2,sprintf('KE = %4.0f meV ',1e3*KE));
set(TT,'fontsize',12)
TT = text(22, -.2,sprintf('PE = %4.1f meV ',1e3*PE));
set(TT,'fontsize',12)
E_tot = KE + PE;
TT = text(22,.2,sprintf( 'E_t_o_t = %4.1f meV ',1e3*(KE+PE)));
set(TT,'fontsize',12)
set(gca,'fontsize',12)
xlabel('nm')
%TT = ylabel( 'f_1','FontName','Symbol','fontsize',12)
TT = ylabel( 'y','FontName','Symbol','fontsize',12)

title('Se2-1')
% The eigenfunction decomposition
angle = zeros(1,40);
for m =1:40
    xeig(m) = m;
    Crl = 0.;
    Cim = 0.;
    for n = 1:NN
        Crl = Crl + prl(n) *phi(n,m);
        Cim = Cim + pim(n) *phi(n,m);
    end
    P_eig(m) = Crl + i *Cim; % The complex value of cn
end

subplot(3,4,3)
bar(xeig,abs(P_eig))
axis( [ 1 20 0 1.1 ])
%TT = text(5,.85,[ '% Eigenfunctions ']);
%set(TT,'fontsize',10)
set(gca,'fontsize',12)
title('Amplitude')
xlabel('Eigenfunction # ')
subplot(3,4,4)
%bar(xeig,angle)
bar(xeig,imag(log(P_eig)))
axis( [ 0 20 -pi pi ])
%TT = text(2,2,[ 'Angle Eigenfunctions ']);
%set(TT,'fontsize',12)
set(gca,'fontsize',12)
title('Phase')
xlabel('Eigenfunction # ')
