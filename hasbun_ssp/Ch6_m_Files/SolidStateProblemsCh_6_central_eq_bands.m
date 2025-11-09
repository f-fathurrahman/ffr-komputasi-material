%copyright by J. E Hasbun and T. Datta
%central_eq_bands.m
%Plots the numeric solutions for the central equation
%for an Nb band system. The eigenvalues are obtained here
%after setting up the matrix of the coefficients.
clear;
clear all;
h=6.62606896e-34;           %Planck's constant (J.s)
hbar=h/2./pi;               %hbar (J.s)
me=9.10938215e-31;          %electron mass (kg)
e=1.602176487e-19;          %electronic charge
ab=1e-10;                   %1 angstom unit of distance
kb=1/ab;                    %wavevector unit
Eb=hbar^2*kb^2/2/me;        %energy unit in joules
Eb_eV=Eb/e;                 %energy unit in eV
%fprintf('Energy unit Eb=%5.3g J, or %5.3f eV\n',Eb,Eb_eV)
a=1.0;            %lattice constant in ab units
g=2*pi/a;
ak=-g/2:0.01:g/2; %k will be chosen parallel to g
ug=2;             %potential magnitude in Eb units
Nb=2;             %number of bands to do
%mg is the value of m in (k-(i-m)g) in the matrix elements
%The mg value chosen seems to work well for symmetry reasons.
if(mod(Nb,2)==0), mg=Nb/2; else mg=(Nb+1)/2; end
aM=zeros(Nb,Nb);  %define the matrix of coefficients
for ik=1:length(ak)
  for i=1:Nb
    %diagonal terms
    if(abs(i-mg) <= 1.e-3)
      aM(i,i)=ak(ik)^2;
    else
      thkg=acos(dot(ak(ik),g)/abs(ak(ik)*g));  %angle between k and g
      aM(i,i)=ak(ik)^2+((i-mg)*g)^2-2*ak(ik)*((i-mg)*g)*cos(thkg);
    end
  end
  for i=1:Nb-1
    aM(i,i+1)=ug;  %complete the tridiagonal matrix
    aM(i+1,i)=ug;
  end
  %get eigenvalues from smallest to highest for each k. Apparently,
  %there is no sorting needed here.
  Eiv(:,ik)=eig(aM);
end
hold on
for i=1:Nb
  plot(ak,Eiv(i,:),'k')
end
xlabel('-g/2 < k < g/2','FontSize',14)
ylabel('\epsilon (\epsilon_b)','FontSize',14)
str=cat(2,num2str(Nb,'%2.0f'),' bands - numeric: a=',...
  num2str(a,'%2.1f'),' Angs, \epsilon_b=',...
  num2str(Eb_eV,'%3.2f'),' eV, g=',num2str(g,'%3.2f'),...
' Angs^{-1}, u_g=',num2str(ug,'%3.1f'),' \epsilon_b');
title(str,'FontSize',12)
axis tight
