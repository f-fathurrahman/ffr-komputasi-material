%copyright by J. E Hasbun and T. Datta
% ch13_metalR.m

% This script plots the reflectance from a metal below, at, and above
% the plasma frequency using the Drude free electron theory for a weakly
% damped system.

% Function definitions
% Equation (13.98)
% The the imaginary part is zero in the limit of weak damping.
% The function is scaled w.r.t to the plasma frequency.

epsr = @(w) (1 - 1/w^2);

% Refractive index, Equation (13.69)

nrefr = @(w) sqrt(epsr(w));

% Reflectance, Equation (13.99)
refl = @(w) (abs((nrefr(w)-1)/(nrefr(w)+1)))^2;

% Plots

figure1 = figure;
axes1 = axes('Parent',figure1,'YMinorTick','on',...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12],...
    'XMinorTick','on');
ylim(axes1,[0 1.2]);
hold(axes1,'all')

hold on;
fplot(@(w)refl(w),[0,3],'LineWidth',2,'-k');

% Create xlabel
xlabel('\omega/\omega_{p}',...
    'FontSize',16,'FontName','Times New Roman');

% Create ylabel
ylabel('Reflectance (R)','FontSize',16,...
    'FontName','Times New Roman');

hold off;
