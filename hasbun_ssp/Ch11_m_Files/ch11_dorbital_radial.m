%copyright by J. E Hasbun and T. Datta
% ch11_dorbital_radial.m

% Hydrogenic radial wavefunction code

% Creating the grid points
rho = linspace(0,20,100);

% Defining the expressions for the radial wavefunction using Table 11.4
% in Chapter 11
% Note: Radial wavefunctions scaled by rho = r/a_o
% Note: Atomic nunmber Z is set to one

Rthreezero = ((2/(3^(3/2))))*(1 - (2/3)*rho + ...
    (2/27)*(rho.*rho)).*exp(-rho/3);               % R_30(r)(3s)
Rthreeone = (4*sqrt(2)/9)*((1/(3^(3/2))))*...
    (1 - rho/6).*rho.*exp(-rho/3);                 % R_31(r)(3p)
Rthreetwo = (4/(27*sqrt(10)))*((1/(3^(3/2))))*rho.*rho...
            .*exp(-rho/3);                         % R_32(r)(3d)

% Radial distribution function. Note the multiplicative r^2 factor in front
% of Rnl(r).

Rthreezerosq = rho.*rho.*Rthreezero.*Rthreezero;          % r^2 R^2_30(r)
Rthreeonesq = rho.*rho.*Rthreeone.*Rthreeone;             % r^2 R^2_31(r)
Rthreetwosq = rho.*rho.*Rthreetwo.*Rthreetwo;             % r^2 R^2_32(r)

% Plotting the n = 3 radial wavefunction
subplot(2,1,1);
plot(rho,Rthreezero,'-k',rho,Rthreeone,'--r',rho,Rthreetwo,'-xb','LineWidth',2);
xlabel('\rho = r/a_{o}')
ylabel('R_{nl}(\rho)')
legend('R_{30}(\rho)','R_{31}(\rho)','R_{32}(\rho)')
title('Hydrogenic Radial Wavefunction, n = 3, l = 0, 1, 2')

% Plotting the radial distrubution functions
subplot(2,1,2)
plot(rho,Rthreezerosq,'o--g',rho,Rthreeonesq,'-r',rho,Rthreetwosq,'-xb',...
    'LineWidth',2);
xlabel('\rho = r/a_{o}')
ylabel('r^2 R^{2}_{nl}(\rho)')
legend('r^2 R_{30}(\rho)','r^2 R_{31}(\rho)','r^2 R_{32}(\rho)');
title('Radial distribution function')
