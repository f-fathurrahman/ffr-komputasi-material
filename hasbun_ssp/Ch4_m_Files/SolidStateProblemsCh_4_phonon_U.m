%copyright by J. E Hasbun and T. Datta
%phonon_U.m
%Calculates the phonon internal energy for a Copper crystal
%based on the Debye Model.
clear;
h=6.62606896e-34;     %Planck'constant (J.s)
hbar=h/2./pi;         %hbar
kB=1.3806504e-23;     %Boltzmann Contant  (J/K)
u=1.660538782e-27;    %atomic mass unit
Na=(1/u)*1.e-3;       %Avogadro's number
Mw=63.546;            %Cu gr/mol. Note 1 gr/mol=1u
Y=0.76e11;            %Cu (cast) Young's mod. in Pa (M. Marder's text)
rho=8890;             %Cu density kg/m^3
vs=sqrt(Y/rho);       %speed of sound
nv=rho/u/Mw;          %N/V
omD=vs*(6*pi^2*nv)^(1/3);  %Debye frequency
thD=hbar*omD/kB;      %Debye temperature
fprintf('N/V=%6.3e 1/m^3, vs=%6.3f m/s\n',nv,vs)
fprintf('Debye freq.=%6.3e 1/m^3, Debye Temp=%6.3f m/s\n',omD,thD)
Uint=@(x) x.^3./(exp(x)-1);              %U integrand
Tol=1.e-3;            %small number for lowest T
T=0:5:thD;            %T variable of U
T(1)=1000*Tol;        %use this value instead of T=0 low limit;
for i=1:length(T)
  xD=thD/T(i);
  %integrate; use small number for lower limit; upper limit is xD
  U(i)=9*Na*kB*T(i)*(T(i)/thD)^3*quad(Uint,Tol,xD);
  UlT(i)=3*pi^4*Na*kB*T(i)^4/thD^3;  %Low T approx
  UhT(i)=3*Na*kB*T(i)-9*Na*kB*thD/8; %High T approx
end
plot(T,U,'k-.','LineWidth',2)
hold on
plot(T,UlT,'kd','MarkerSize',4)
plot(T,UhT,'ko','MarkerSize',4)
legend('Numeric','Low T Approx','High T Approx',4);
axis([0 max(T) -100 max(U)])
xlabel('T (Kelvin)')
ylabel('U (J/mol)')
