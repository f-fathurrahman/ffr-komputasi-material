%copyright by J. E Hasbun and T. Datta
% ch13_nwkappawmetal.m

% This script plots the refractive index (n), extinction coefficient
% kappa, and the reflectance derived within the Drude theory in the case
% when damping is not negligible. The plasma frequency and damping
% parameters are for Na (metal). Data sourced from Reference [13.3],
% Table 11.2.

% Constants
omp = 5.914;   % Na plasma frequency (eV)
hgamma = 0.0198512; % Na damping frequency (eV)

% Function definitions
% Frequency dependent real part of the dielectric function.
% Compare with Equation (13.96).
realpart = @(w) 1 - omp^2/(hgamma^2 + w^2);

% Frequency dependent imagninary part of the dielectric function.
% Compare with Equation (13.97).
impart = @(w) (omp^2)*hgamma/((w)*(w^2 + hgamma^2));

% Refractive index (n) and extinction coefficient (\kappa)
nw = @(w) sqrt(0.5*(sqrt((1+realpart(w))^2 + impart(w)^2) + 1 ...
    + realpart(w)));
kappaw = @(w) sqrt(0.5*(sqrt((1+realpart(w))^2 + impart(w)^2)...
    - 1 - realpart(w)));

% Reflectance, Equation (13.99)

refl = @(w) ((1 - nw(w))^2 + kappaw(w)^2)/((1 + nw(w))^2 + kappaw(w)^2);

% Plots

figure1 = figure;
axes1 = axes('Parent',figure1,'YMinorTick','on',...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12],...
    'XMinorTick','on');
ylim([0,2]);
hold(axes1,'all')

hold on;

fplot(@(w)nw(w),[0,10],'LineWidth',2,'--k');
fplot(@(w)kappaw(w),[0,10],'LineWidth',2,'-or');
fplot(@(w)refl(w),[0,10],'LineWidth',2,'-b');

% Create xlabel
xlabel('Energy,(eV)',...
    'FontSize',16,'FontName','Times New Roman');

% Create ylabel
ylabel('Optical Constants and Reflectance', 'FontSize',16,...
    'FontName','Times New Roman');

% % Create legend
legend1 = legend('n(\omega)','\kappa(\omega)','R(\omega)');
set(legend1,'FontSize',14,'FontName','Times New Roman','show');

hold off;
