%copyright by J. E Hasbun and T. Datta
%orbitals_vs_energy.m
%The number of orbitals per energy versus energy
%as the product of the 3-Dim dos times the Fermi-Dirac
%distribution function is plotted versus energy.
clear
h=6.62606896e-34;           %Planck'constant (J.s)
hbar=h/2./pi;               %hbar
me=9.10938215e-31;          %electron mass (kg)
e=1.602176487e-19;          %electronic charge
Eb=e;                       %energy in Joules (=1eV)
kB=1.3806504e-23;           %Boltzmann Constant (J/K)
Tb=Eb/kB;                   %temperature unit
Lb=(hbar^2/2/me/Eb)^(1/2);  %length unit
fprintf('Energy unit: %3.2f eV,\n',Eb/e)
fprintf('Temperature unit: %4.3g K\n',Tb)
fprintf('Length unit: %4.3g m\n',Lb)
%temperature variable in units of Tb
Tv=[0.01,0.5];
Emax=4;
Emin=0;
Es=(Emax-Emin)/500;
mu=2.5;
eps=0:Es:Emax;
fe= @(eps,mu,T) 1./(exp((eps-mu)./T)+1);  %FD distribution definition
De= @(eps) eps.^(1/2)/2/pi^2;             %Density of states 3D definition
for j=1:length(Tv)
  T=Tv(j);
  subplot(1,2,j)
  hold on
  for i=1:length(eps)
    fd(i)=fe(eps(i),mu,T);
    dos(i)=De(eps(i));
    ne(i)=fd(i)*dos(i);                  %orbitals/energy
  end
  plot(eps,fd,'k-.')
  plot(eps,2*pi^2*dos,'k:','LineWidth',2)
  plot(eps,2*pi^2*ne,'k-')
  legend('f(\epsilon)','2\pi^2D(\epsilon)',...
    '2\pi^2f(\epsilon)D(\epsilon)',0)
  xlabel('\epsilon (eV)','Fontsize',14)
  str1='f(\epsilon), 2\pi^2D(\epsilon), 2\pi^2f(\epsilon)D(\epsilon)';
  ylabel(str1,'Fontsize',14)
  str2=cat(2,'\mu=',num2str(mu,'%3.1f'),'\epsilon_b, T=',...
    num2str(T,'%3.1f'),' T_b');
  title(str2)
end
