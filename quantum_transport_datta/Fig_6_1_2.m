clear all
close all

z = zeros(5);
Z = zeros(10);

%Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
a = 2.45e-10*4/sqrt(3);
m = 9.110e-31;

d1 = [1 1 1]/4;
d2 = [1 -1 -1]/4;
d3 = [-1 1 -1]/4;
d4 = [-1 -1 1]/4;

%sp3s* model parameters
soa = 0.3787/3;
soc = 0.0129/3;
Esa = -8.3431;
Epa = 1.0414;
Esc = -2.6569;
Epc = 3.6686;
Esea = 8.5914;
Esec = 6.7386;
Vss = -6.4513;
Vpasc = -5.7839;
Vpasec = -4.8077;
Vsapc = 4.4800;
Vseapc = 4.8422;
Vxx = 1.9546;
Vxy = 5.0779;

%Conduction band effective mass model parameters
Ec = 1.55;
meff = 0.12*m;
Nt = 101;
kk = 1*linspace(0,1,Nt);

l=0.5; m=0.5; n=0.5; %L-direction
%l=1;m=0;n=0;%X-direction

for Nk = 1:Nt
    k = 2*pi*kk(Nk)*[l m n];
    %sp3s* model
    p1 = exp(i*sum(k.*d1));
    p2 = exp(i*sum(k.*d2));
    p3 = exp(i*sum(k.*d3));
    p4 = exp(i*sum(k.*d4));
    g0 = (p1+p2+p3+p4)/4;
    g1 = (p1+p2-p3-p4)/4;
    g2 = (p1-p2+p3-p4)/4;
    g3 = (p1-p2-p3+p4)/4;
    a1 = diag([Esa Epa Epa Epa Esea]);
    A1 = [a1 z;
          z a1];
    a2 = diag([Esc Epc Epc Epc Esec]);
    A2 = [a2 z;
          z a2];
    b = [Vss*g0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0;
         Vpasc*g1 Vxx*g0 Vxy*g3 Vxy*g2 Vpasec*g1;
         Vpasc*g2 Vxy*g3 Vxx*g0 Vxy*g1 Vpasec*g2;
         Vpasc*g3 Vxy*g2 Vxy*g1 Vxx*g0 Vpasec*g3;
         0 Vseapc*conj(g1) Vseapc*conj(g2) Vseapc*conj(g3) 0];
    B = [b z;
         z b];
    h = [a1 b;
         b' a2];
    H = [A1 B;
         B' A2];
    aso = soa*[0 0 0 0 0 0 0 0 0 0;
               0 0 -i 0 0 0 0 0 1 0;
               0 i 0 0 0 0 0 0 -i 0;
               0 0 0 0 0 0 -1 i 0 0;
               0 0 0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0 0 0;
               0 0 0 -1 0 0 0 i 0 0;
               0 0 0 -i 0 0 -i 0 0 0;
               0 1 i 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0 0 0];
    cso=soc*[0 0 0 0 0 0 0 0 0 0;
             0 0 -i 0 0 0 0 0 1 0;
             0 i 0 0 0 0 0 0 -i 0;
             0 0 0 0 0 0 -1 i 0 0;
             0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0;
             0 0 0 -1 0 0 0 i 0 0;
             0 0 0 -i 0 0 -i 0 0 0;
             0 1 i 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0];
    H = H + [aso Z; Z cso];
    %
    [V,D]=eig(H);
    eiglst = sum(D);
    E(Nk,:) = sort(real(eiglst));
    %Conduction band effective mass model
    Em(Nk) = Ec + ((hbar^2)*sum(k.*k)/(2*meff*q*(a^2)));
end

kk = -kk; %L-direction
hold on
h1=plot(kk,E,'b');
h2=plot(kk,Em,'b--');
axis([-1 1 -3 3])
set(h1,'linewidth',[1.0])
set(h2,'linewidth',[2.0])
set(gca,'Fontsize',[24])
xlabel('ka (fraction of maximum value')
ylabel('Energy (eV)')
grid on
