%copyright by J. E Hasbun and T. Datta
%alloy_CuNi_cond.m
%Case of 293 K
%To use the effective medium conductivity formula to obtain
%the conductivity versus concentration.
%We work with the Cu(x)Ni(1-x) alloy.
%Cu(x)Ni(1-x) - g in units of 10^8/(ohm-meter)
clear, clc
%recommended data from C. Y. Ho et al.,
%J. Phys. Chem. Ref. Data V12 (2) 183 (1983)
x1d=[0.000,0.005,0.010,0.030,0.050,0.100,0.150,0.200,...
  0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,...
  0.700,0.750,0.800,0.850,0.900,0.950,0.970,0.990,0.995,1.000];
g1d=[0.144,0.133,0.124,0.097,0.080,0.056,0.043,0.035,0.028,...
  0.024,0.022,0.021,0.020,0.020,0.021,0.022,0.024,0.027,...
  0.032,0.039,0.050,0.072,0.129,0.188,0.351,0.448,0.596];
Nd=length(x1d);
plot(x1d,g1d,'ko')
hold on
g1=g1d(Nd);    %pure copper g
g2=g1d(1);     %pure nickel g
x=0:0.05:1;
Nx=length(x);
for i=1:Nx
  y=1-x(i);
  gbar(i)=x(i)*g1+y*g2; %VCA
end
plot(x,gbar,'k--') %VCA Cu-Ni
