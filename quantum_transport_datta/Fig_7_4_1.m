clear all
close all

z = zeros(5);
Z = zeros(10);

% Constants (all MKS, except energy which is in eV)
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
Epc=3.6686;Esea=8.5914;Esec=6.7386;
Vss=-6.4513;Vpasc=-5.7839;Vpasec=-4.8077;
Vsapc=4.4800;Vseapc=4.8422;Vxx=1.9546;Vxy=5.0779;

%Valence band Luttinger-Kohn parameters
Ev=-.1;del=.3;g1=6.85;g2=2.1;g3=2.9;
t1=(hbar^2)*g1/(2*m*q*(a^2));
t2=(hbar^2)*g2/(2*m*q*(a^2));
t3=(hbar^2)*g3/(2*m*q*(a^2));

Nt=101;kk=1*linspace(0,1,Nt);
l=1;m=0;n=0;%X-direction
l=0.5;m=0.5;n=0.5;%L-direction

for Nk=1:Nt
    k=2*pi*kk(Nk)*[l m n];
    %sp3s* model
    p1 = exp(i*sum(k.*d1)); p2 = exp(i*sum(k.*d2));
    p3 = exp(i*sum(k.*d3)); p4 = exp(i*sum(k.*d4));
    g0 = (p1 + p2 + p3 + p4)/4; g1 = (p1 + p2 - p3 - p4)/4;
    g2 = (p1 - p2 + p3 - p4)/4; g3 = (p1 - p2 - p3 + p4)/4;
    a1 = diag([Esa Epa Epa Epa Esea]);A1=[a1 z;z a1];
    a2 = diag([Esc Epc Epc Epc Esec]);A2=[a2 z;z a2];
    b = [Vss*g0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0;
         Vpasc*g1    Vxx*g0 Vxy*g3 Vxy*g2 Vpasec*g1;
         Vpasc*g2 Vxy*g3 Vxx*g0 Vxy*g1 Vpasec*g2;
         Vpasc*g3 Vxy*g2 Vxy*g1 Vxx*g0 Vpasec*g3;
         0 Vseapc*conj(g1) Vseapc*conj(g2) Vseapc*conj(g3) 0];
    B = [b z;z b];
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
    cso = soc*[0 0 0 0 0 0 0 0 0 0;
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
    [V,D] = eig(H);
    eiglst = sum(D);
    E(Nk,:) = sort(real(eiglst));

    %Valence band Luttinger-Kohn model
    P = Ev + (t1*sum(k.*k));
    Q = t2*((k(1)^2)+(k(2)^2)-(2*(k(3)^2)));
    R = -(sqrt(3)*t2*((k(1)^2)-(k(2)^2)))+(i*2*t3*sqrt(3)*k(1)*k(2));
    S = 2*t3*sqrt(3)*((k(1)-(i*k(2)))*k(3));

    H4 = -[P+Q -S R 0;
          -S' P-Q 0 R;
           R' 0 P-Q S;
           0 R' S' P+Q];
    [V,D] = eig(H4);
    eiglst = sum(D);
    ELK4(Nk,:) = sort(real(eiglst));

    H6 = -[P+Q -S R 0 -S/sqrt(2) sqrt(2)*R;
          -S' P-Q 0 R -sqrt(2)*Q sqrt(1.5)*S;
           R' 0 P-Q S sqrt(1.5)*S' sqrt(2)*Q;
           0 R' S' P+Q -sqrt(2)*R' -S'/sqrt(2);
          -S'/sqrt(2) -sqrt(2)*Q' sqrt(1.5)*S -sqrt(2)*R  P+del 0;
          sqrt(2)*R' sqrt(1.5)*S' sqrt(2)*Q' -S/sqrt(2) 0 P+del];
    [V,D] = eig(H6);
    eiglst = sum(D);
    ELK6(Nk,:) = sort(real(eiglst));
end

kk = -kk; %L-direction

hold on
h1=plot(kk,E,'b');
%h2=plot(kk,ELK4,'b--');% Fig.6.4.1
h2=plot(kk,ELK6,'b--');% Fig.6.4.2
set(h1,'linewidth',[2.0])
set(h2,'linewidth',[3.0])
set(gca,'Fontsize',[24])
xlabel(' ka (fraction of maximum value)')
ylabel(' Energy (eV)')
axis([-1 1 -2 3])
grid on