%copyright by J. E Hasbun and T. Datta
% ch13_drudelorentz.m

% This script plots the real and imaginary parts of the susceptibility
% derived within the Drude-Lorentz theory.

% Constants
ompscl = 2;   % scaled plasma frequency, \omega_p/\omega_o
gammascl = 0.25; % scaled damping rate, \gamma/\omega_o

% Function definitions
% Frequency dependent real part of the susceptibility (scaled).
% Compare with Equation (13.125). The frequency w is scaled w.r.t omega_{o}

realpart = @(w) (ompscl^2)*(1 - w^2)/((1 - w^2)^2+gammascl^2*w^2);

% Frequency dependent imaginary part of the susceptibility (scaled).
% Compare with Equation (13.126).
impart = @(w) (ompscl^2)*(gammascl*w)/((1 - w^2)^2+gammascl^2*w^2);

% Plots

figure1 = figure;
axes1 = axes('Parent',figure1,'YMinorTick','on',...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12],...
    'XMinorTick','on');
hold(axes1,'all')

hold on;
fplot(@(w)realpart(w),[0,3],'LineWidth',2,'-k');
fplot(@(w)impart(w),[0,3],'LineWidth',2,'--r');

% Create xlabel
xlabel('Scaled frequency, \omega/\omega_{o}',...
    'FontSize',14,'FontName','Times New Roman');

% Create ylabel
ylabel('Susceptibility, \epsilon(\omega)','FontSize',14,...
    'FontName','Times New Roman');

% Create legend
legend('\Ree[\epsilon(\omega)] - 1 \equiv \epsilon_{1}(\omega)',...
    '\Imm[\epsilon(\omega)]\equiv \epsilon_{2}(\omega)');
legend(axes1,'show');

% Create textbox
annotation(figure1,'textbox',...
    [0.378380952380953 0.602947368421053 0.0489999999999994 0.118],...
    'String',{'\gamma'},...
    'FontSize',30,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create doublearrow
annotation(figure1,'doublearrow',[0.353571428571429 0.417857142857143],...
    [0.632 0.634]);

% Create line
annotation(figure1,'line',[0.354761904761905 0.361904761904762],...
    [0.924 0.108],'LineStyle',':','LineWidth',2);

% Create line
annotation(figure1,'line',[0.388095238095238 0.383333333333333],...
    [0.107270676691729 0.926],'LineStyle','--','LineWidth',2);

% Create line
annotation(figure1,'line',[0.416666666666667 0.415476190476191],...
    [0.11 0.924],'LineStyle',':','LineWidth',2);

% Create line
annotation(figure1,'line',[0.131897711978466 0.909825033647376],...
    [0.377594249201278 0.381789137380192],'LineWidth',2);

% Create line
annotation(figure1,'line',[0.129761904761905 0.907142857142857],...
    [0.409 0.414],'LineStyle','--','LineWidth',2);

hold off;
