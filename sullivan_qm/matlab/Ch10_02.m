% Se2d_pert.m.  Adds a perturbation gradually using a Hanning window.

clear all

NN = 50;              % Number of points in the problem space.
MM = 50;
hbar = 1.054e-34;     % Plank's constant
melec = 9.1e-31;      % Mass of an electron
eV2J = 1.6e-19;       % Energy conversion factors
J2eV = 1/eV2J;

del_x = 2e-10;        % The cells size
dt = .05e-15;           % Time steps 
ra = (0.5*hbar/melec)*(dt/del_x^2) % ra must be < .1
DX = del_x*1e9;       % Cell size in nm.
XX = (DX:DX:DX*NN);   % Length in nm for plotting
YY = (DX:DX:DX*MM);   % Length in nm for plotting

NC = NN/2;            % Starting position of the pulse
MC = NN/2;

%  --- Hanning Window

han = zeros(1,3e4);
del_T = dt*1e15;
Time = (0:del_T:del_T*(3e4-1));
for n=1:10000
   han(n) = 0.5*(1-cos(2*pi*n/9999)) ;
   %han(n) = exp(-0.5*((n-5000)/1500)^2);
end

% --- Add the purturbing potential

V = zeros(NN,MM);
Vp = input('Vpert (eV)  -->');

for n=2:NC
   for m=2:MC
      V(n,m) = Vp;
   end 
end

subplot(2,2,3)
mesh(YY,XX,1e3*V)
view(30,30)
axis( [ 0 DX*MM 0 DX*NN 0. 1.2e3*Vp ])
xlabel('Y (nm)')
ylabel('X (nm)')
zlabel('H_p (meV)')
title('Se2d-pert')
set(gca,'fontsize',12)

subplot(3,3,1)
plot(Time,han,'k')
axis( [ 0 500 0 1 ])
xlabel('fs')
ylabel('Han(t)')
set(gca,'fontsize',12)

%---- Input -----

prl = zeros(NN,MM);
pim = zeros(NN,MM);
ptot = 0.;
for n=2:NN-1
   for m=2:MM-1
      %prl(n,m) = sin(pi*n/NN)*sin(2*pi*m/MM);
      prl(n,m) =  sin(pi*n/NN)*sin(2*pi*m/MM) - sin(2*pi*n/NN)*sin(1*pi*m/MM);  % The Good state
      pim(n,m) = 0.;
      ptot = ptot + pim(n,m)^2 + prl(n,m)^2;
   end 
end
ptot

ptot1 = 0.;
for n=1:NN
   for m=1:MM
      prl(n,m) = prl(n,m)/sqrt(ptot);
      pim(n,m) = pim(n,m)/sqrt(ptot);
      ptot1 = ptot1 + prl(n,m)^2 + pim(n,m)^2;
   end 
end
ptot1

%saveas(gcf,'pert.png')      % This saves the picture to a file

T = 0;
n_step = 1;

%while n_step > 0
n_step = input('How many time steps  -->');

% -----------This is the core FDTD program -------------
for iT=1:n_step

   T = T + 1;

   for m=2:MM-1
      for n=2:NN-1
         prl(n,m) = prl(n,m) - ra*(-4*pim(n,m) + pim(n,m-1) + pim(n,m+1)  ...
                  + pim(n-1,m) + pim(n+1,m) ) ...
                  + han(T)*(dt/hbar)*eV2J*V(n,m)*pim(n,m);
      end 
   end 

   for m=2:MM-1
      for n=2:NN-1
         pim(n,m) = pim(n,m) + ra*(-4*prl(n,m) + prl(n,m-1)  + prl(n,m+1)  ...
                    + prl(n-1,m)  + prl(n+1,m) ) ... 
                    - han(T)*(dt/hbar)*eV2J*V(n,m)*prl(n,m);
      end
   end

end
T
han(T)
% ------------------------

% Check normalization
ptot = 0.;
for m=1:MM
   for n=1:NN
      psi(n,m) = prl(n,m) + i*pim(n,m);
      ptot = ptot +  prl(n,m)^2 + pim(n,m)^2;
   end
end
ptot 

% Calculate the expected values
ke = 0. + j*0.;
pe = 0.;
for m=2:MM-1
   for n=2:NN-1
      lap_p = psi(n+1,m) - 4*psi(n,m) + psi(n-1,m) + psi(n,m-1) + psi(n,m+1);
      ke = ke + lap_p*psi(n,m)';
      pe = pe + han(T)*V(n,m)*psi(n,m)*psi(n,m)';
   end
end
KE = -J2eV*((hbar/del_x)^2/(2*melec))*real(ke);

subplot(2,2,2)
mesh(YY,XX,prl) 
view(30,30)
axis( [ 0 DX*MM 0 DX*NN -.03 .03 ])
% TT = text(0,0,.05,sprintf('%4.0f fs',T*dt*1e15));
% set(TT,'fontsize',12)
%TT = text(0,0,.07,sprintf('KE = %5.2f meV  PE = %5.2f meV',1e3*KE,1e3*pe));
%set(TT,'fontsize',12)
TT = zlabel('Y')
set(TT,'FontName','symbol')
set(TT,'FontName','symbol')
set(gca,'fontsize',12)
xlabel('Y (nm)')
ylabel('X (nm)')
title('Se2d-pert')

subplot(3,3,9)
contour(YY,XX,prl)
TT = text(1.5,8,sprintf('%4.0f fs',T*dt*1e15));
set(TT,'fontsize',12)
xlabel('Y (nm)')
ylabel('X (nm)')
set(gca,'fontsize',11)

%saveas(gcf,'pert.png')      % This saves the picture to a file

%end % while nstep

