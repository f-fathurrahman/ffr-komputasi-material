%Example for the calculation of a photon's wavelength versus
%photon energy: a log-log plot
clear;
hc=12398;              %h*c in eV-angstroms
N=100;                 %steps to use in the plot
%Ek = particles' relativistic kinetic energy
Eki=5e2; Ekf=1e5;      %photon E range (eV): initial, final
Eks=(Ekf-Eki)/(N-1);   %step size
Ek=Eki:Eks:Ekf;
lambda=hc./Ek;
loglog(Ek,lambda,'k')  %log-log plot
axis tight
legend('Photons (E_k in eV)',0);
xlabel('E_k')
ylabel('\lambda (Å)')
