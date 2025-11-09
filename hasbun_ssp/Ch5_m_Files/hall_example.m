%copyright by J. E Hasbun and T. Datta
%hall_example.m
%Simple illustration of the classical transverse and
%logitudinal resistivities along with their quantum mechanical
%versions.
clear
e=1.602176487e-19;        %electronic charge
h=6.62606896e-34;         %Planck'constant (J.s)
fxQ=h/e;                  %quantized unit of flux (webers)
rhoQ=fxQ/e;               %quantized resistivity (ohm) (in 2D)
A=1.e-4;                  %sample area m^2
N=3.e10;                  %number of electrons
Bmin=0; Bs=0.25; Bmax=10; %range of B
B=Bmin:Bs:Bmax;
nB=length(B);
for i=1:nB
  rhoH(i)=B(i)*A/N/e;  %classical rho_H in Ohm (classical)
end
rhoL=0.015*rhoH(nB);   %class. rho_L in Ohm - let it be some constant
%Quantized Hall resistance
%Need nu (the Landau level fill factor), use N=nu*e*B*A/h=nu*B*A/fxQ
nu_min=N*fxQ/Bmax/A;    %minimum value of nu
nu_max=N*fxQ/Bs/A;      %maximum value of nu
nu_low=ceil(nu_min);    %round the value
nu_high=ceil(nu_max);   %highest nu for rho to reach rhoH(nB)
nu=nu_low:1:nu_high;    %vary nu in integer values
nNu=length(nu);
for j=1:nNu
  rhoHQ(j)=rhoQ/nu(j);  %quantum rho_H=(h/e^2)/nu
  %Find the B field values of rho_L=0, use def: BA=N*(h/e)/nu
  BQ(j)=N*fxQ/A/nu(j);
end
plot(B,rhoH/1e4,'k')    %classical rho_H
hold on
%classical rho_L
line([B(1) B(nB)],[rhoL/1e4 rhoL/1e4],'Color','k','LineStyle',':',...
  'LineWidth',2)
xlabel('B_z (Tesla)')
ylabel('\rho (10^4 \Omega)')
plot(BQ,rhoHQ/1e4,'ksq') %quantum rho_H
plot(BQ,0.0,'ko')        %quantum rho_L zero positions
for j=1:nNu
  str=cat(2,'  \nu=',num2str(nu(j),'%i'));
  text(BQ(j),rhoHQ(j)/1e4,str)
end
axis([0 6*BQ(nNu) 0 6*rhoHQ(nNu)/1e4])
legend('Classical \rho_H','Classical \rho_L','Quantum \rho_H plateau positions',...
  'Quantum \rho_L zero positions',0)
str=cat(2,'Classical, Integer Quantum Hall Effect Example: ',...
    'N_s=',num2str(N/(A*1e4),'%6.2e'),'cm^{-2}');
title(str)
