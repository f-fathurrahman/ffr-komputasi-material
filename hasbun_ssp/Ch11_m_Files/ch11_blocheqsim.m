%copyright by J. E Hasbun and T. Datta
% ch11_blocheqsim.m

% This script simulates the Bloch equations in the presence of
% logitudinal and transverse relaxation

% Solving the Bloch equations
% dMx_dt(1) = write all components
% dMy_dt(2) = write all components
% dMz_dt(3) = write all components

%initial conditions
Xo = [1;0;0];
%timespan
tspan = [0,5];
% defining simulation parameters
omega= 50;tone=1;ttwo=1;mzero=1;
%set an error for solving coupled ODE
options=odeset('RelTol',1e-6);
magparam = @(t,x)mag(t,x,omega,tone,ttwo,mzero);
%call the solver
[t,X] = ode45(magparam,tspan,Xo,options);
%plot the results
figure(1)
hold on;
plot(t,X(:,1),'r');plot(t,X(:,2),':');plot(t,X(:,3),'k-');
legend('M_x(t)','M_y(t)','M_z(t)');ylabel('x');xlabel('t')
figure(2)
plot(X(:,1),X(:,2));

% Create a separate function file ch11_mag.m

function dm = mag(~,m,omega,tone,ttwo,mzero)
%a function which returns a rate of change vector
dm = zeros(3,1);
dm(1)=  omega*m(2)- m(1)/ttwo;
dm(2)= -omega*m(1)- m(2)/ttwo;
dm(3)=  (mzero - m(3))/tone;
end
