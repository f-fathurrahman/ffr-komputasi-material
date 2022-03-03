% 1d FDTD simulation of psi

clear all; close all

% Number of points in the problem space.
NN = 400;

% Planck constant
hbar = 1.054e-34;

% Free space mass of an electron
m0 = 9.1e-31;

% Effective mass: Si is 1.08, Ge is 0.067, GaAs is 0.55
meff = 1.0;

% Mass of an electron
melec = meff*m0;

% Charge of an electron
ecoul = 1.6e-19;

% Dielectric of free space
epsz = 8.85e-9;

% Energy conversion factors
eV2J = 1.6e-19;
J2eV = 1/eV2J;

% The cell size
del_x = 0.1e-9;

% Time steps
dt = 2e-17;

% ra must be < 0.15
ra = (0.5*hbar/melec)*(dt/del_x^2);

% Cell size in nm.
DX = del_x*1e9;

% Length in nm for plotting
XX = (DX:DX:DX *NN);

% Specify the potential
V=zeros(1,NN);

% Barrier
for n=NN/2:NN/2+50
    V(n) = .15*eV2J;
end

% Semiconductor conduction band
%for n=1:NN/2
%    %V(n) = .1*eV2J;
%end

%for n=NN/2+1:NN
%    %V(n) = .2*eV2J;
%end

% Electric ﬁeld
%for n=1:NN
%    V(n) = -(0.2/400)*(n-400)*eV2J;
%end

% Initialize a sine wave in a gaussian envelope

% Pulse wavelength
lambda = 50;
%lambda = 25;

sigma = 50; % Pulse width

nc = 150; % Starting position

prl = zeros(1,NN); % The real part of the state variable
pim = zeros(1,NN); % The imaginary part of the state variable
ptot = 0.0;
for n=2:NN-1
    prl(n) = exp(-1.*((n-nc)/sigma)^2)*cos(2*pi*(n-nc)/lambda) ;
    pim(n) = exp(-1.*((n-nc)/sigma)^2)*sin(2*pi*(n-nc)/lambda) ;
    ptot = ptot + prl(n)^2 + pim(n)^2;
end
pnorm = sqrt(ptot); % Normalization constant

% Normalize and check
ptot = 0.;
for n=1:NN
    prl(n) = prl(n)/pnorm;
    pim(n) = pim(n)/pnorm;
    ptot = ptot + prl(n)^2 + pim(n)^2;
end

fprintf('ptot = %f\n', ptot) % This should have the value 1

T = 0;
n_step = 500;
%n_step = input('How many time steps -->');

% -----------This is the core FDTD program ------------
for m=1:n_step
    T = T + 1;
    for n=2:NN-1
        prl(n) = prl(n) - ra*(pim(n-1) -2*pim(n) + pim(n+1)) + (dt/hbar)*V(n)*pim(n);
    end
    for n=2:NN-1
        pim(n) = pim(n) + ra*(prl(n-1) -2*prl(n) + prl(n+1)) - (dt/hbar)*V(n)*prl(n);
    end
end

% ------------------------------------------------
% Calculate the expected values
PE = 0.0;
for n=1:NN
    % Write as a complex function
    psi(n) = prl(n) + i*pim(n);
    PE = PE + psi(n)*psi(n)'*V(n);
end
    
% This checks normalization
fprintf('Check norm: %f\n', psi*psi')
    
PE = PE*J2eV; % Potential energy

ke = 0.0 + j*0.;
for n=2:NN-1
    lap_p = psi(n+1) - 2*psi(n) + psi(n-1);
    ke = ke + lap_p*psi(n)';
end
% Kinetic energy
KE = -J2eV*((hbar/del_x)^2/(2*melec))*real(ke);

%subplot(2,1,1)
plot(XX,prl,'k')
hold on
plot(XX,pim,'-.k')
plot(XX,J2eV*V,'--k')
hold off
axis( [ 1 DX*NN -0.2 0.3 ])
TT = text(5, 0.15, sprintf('%7.0f fs',T*dt*1e15));
set(TT,'fontsize',12)
TT = text(5, -0.15, sprintf('KE = %5.3f eV',KE));
set(TT,'fontsize',12)
TT = text(25, -0.15, sprintf('PE = %5.3f eV',PE));
set(TT,'fontsize',12)
TT = text(25, 0.13, sprintf('E_t_o_t = %5.3f eV',KE+PE));
set(TT,'fontsize', 12)
xlabel('nm')
set(gca,'fontsize',12)
title('Se1-1')