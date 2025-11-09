%copyright by J. E Hasbun and T. Datta
%mat_sigmaVs_conc.m
%plots the data of some materials' conductivities vs their carrier
%concentrations
clear
%Using typical conductivity values and intrinsic carrier concentrations
mat={'GaAs' 'Si' 'Ge' 'InSb' 'Sb' 'As' 'Na' 'Cu'};
sig=[2.55e-07 4.59e-04 2.23 2.00e+04...
   2.56e6 3.85e6 2.11e+07 5.88e+07];
nv=[1.79e+12 1.45e+16 2.4e+19 1.60674e+22...
  7e25 3e26 2.60e+28 8.45e+28];
%
loglog(nv,sig,'ko','MarkerSize',5)
hold on
nm=length(mat);
text(nv*(1-0.99),sig*2,mat,'FontSize',10)
line([nv(1) nv(end)],[1e5 1e5],'Color','k',...
  'LineStyle','--','LineWidth',2);   %metals above this line
line([nv(1) nv(end)],[1e-6 1e-6],'Color','k',...
  'LineStyle',':','LineWidth',2);    %insulators below this line
xlabel('Carrier Concentration (m^{-3})','FontSize',14)
ylabel('Conductivity ( \Omega\cdot m)^{-1}','FontSize',14)
