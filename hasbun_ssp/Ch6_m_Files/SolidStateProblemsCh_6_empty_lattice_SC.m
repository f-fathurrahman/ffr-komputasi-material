%copyright by J. E Hasbun and T. Datta
%empty_lattice_SC.m
%The empty lattice approximation is used for the
%simple cubic system in the reduced zone scheme
clear
%SC
h=6.62606896e-34;           %Planck'constant (J.s)
hbar=h/2./pi;               %hbar (J.s)
me=9.10938215e-31;          %electron mass (kg)
e=1.602176487e-19;          %electronic charge
ab=1e-10;                   %1 angstom unit of distance
kb=1/ab;                    %wavevector unit
Eb=hbar^2*kb^2/2/me;        %energy unit in joules
Eb_eV=Eb/e;                 %energy unit in eV
%fprintf('Energy unit Eb=%5.3g J, or %5.3f eV\n',Eb,Eb_eV)
a = 2;           %Lattice constant in ab units,
a1=a*[1,0,0];    %using the cubic cell vectors
a2=a*[0,1,0];
a3=a*[0,0,1];
%Reciprocal lattice vectors follow
Vt=dot(a1,cross(a2,a3));  %system's unit cell volume
b1=2*pi*cross(a2,a3)/Vt;
b2=2*pi*cross(a3,a1)/Vt;
b3=2*pi*cross(a1,a2)/Vt;
vmax=1;                      %max v integer (determines bands)
dir=[1,0,0];                 %BZ direction chosen
%max k's on the first BZ - Simple Cubic
Gs=b1+b2+b3;                 %make a small G
kx_max=Gs(1)/2;              %max kx from Gs bisector
ky_max=Gs(2)/2;              %max ky  "    "
kz_max=Gs(3)/2;              %max kz  "    "
kmax=[kx_max,ky_max,kz_max];
%chosen k direction
dirk=[dir(1)*kmax(1),dir(2)*kmax(2),dir(3)*kmax(3)];
x=-1:0.1:1;                  %direction step size
hold on
for v1=-vmax:vmax
  for v2=-vmax:vmax
    for v3=-vmax:vmax
      G=[v1*b1+v2*b2+v3*b3]; %reciprocal lattice vector
      for i=1:length(x)
        kdir=x(i)*dirk;      %step in the chosen direction
        eps(i)=((kdir(1)+G(1))^2+(kdir(2)+G(2))^2+(kdir(3)+G(3))^2);
      end
      plot(x,eps,'k')
    end
  end
end
epm=(Gs(1)^2+Gs(2)^2+Gs(3)^2); %for plotting purposes
axis([-1 1 0 0.75*epm])
xlabel('k_x','FontSize',14)
ylabel('\epsilon (\epsilon_b)','FontSize',14)
str=cat(2,'Reduced BZ Scheme, a=',num2str(a,'%3.2f'),' Angs, ',...
  '\epsilon_b=',num2str(Eb_eV,'%3.2f'),' eV');
title(str,'FontSize',13)
set(gca,'XTick',[-1,1])
set(gca,'XTickLabel','')
set(gca,'XTickLabel',{'-pi/a','pi/a'})
