clear all
close all

t = 3;
m = 800; % Use 200 and 800 for two plots
a0 = 0.14;
D = 2*m*a0*sqrt(3)/(2*pi);
Eg = 2*t*0.14/D;
c = pi*D;
L = 1;
nu0 = round(2*m/3);a=3*a0/2;
E = linspace(0,0.25,101);
DG = (2*c*L/(2*pi*a*a*t*t))*E;
DN = zeros(1,101);
for nu = nu0-100:nu0+100
    Ek = ((t*2*pi/sqrt(3))*((3*nu/(2*m))-1))+(i*1e-12);
    DN = DN + ((2*L/(pi*a*t))*abs(real(E./(sqrt((E.^2)-(Ek^2))))));
end

hold on
h1=plot(DG,E,'bx');
h2=plot(DN,E,'b');
hold on
axis([0 50 0 0.25]);
set(h1,'linewidth',[1.0])
set(h2,'linewidth',[2.0])
set(gca,'Fontsize',[24])
xlabel('D(E) (per eV per nm)')
ylabel('Energy (eV)')
grid on
