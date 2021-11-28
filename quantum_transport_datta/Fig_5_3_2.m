clear all
close all

Esa = -8.3431;
Epa = 1.0414;
Esc = -2.6569;
Epc = 3.6686;
Esea = 8.5914;
Esec = 6.7386;
Vss = -6.4513;
Vxx = 1.9546;
Vxy = 5.0779;
Vsapc = 4.4800;
Vpasc = 5.7839;
Vseapc = 4.8422;
Vpasec = 4.8077;

%Either of the following choices for d1,d2,d3 and d4 should give the same result.
d1=[1 1 1]/4;
d2=[1 -1 -1]/4;
d3=[-1 1 -1]/4;
d4=[-1 -1 1]/4;
d1=[0 0 0]/2;
d2=[0 -1 -1]/2;
d3=[-1 0 -1]/2;
d4=[-1 -1 0]/2;

l = 1; m = 1; n = 1; kmax = pi; Nt = 21; %L-direction

%l=1;m=0;n=0;kmax=2*pi;Nt=21;%X-direction

for Nk = 1:Nt
    k = [l m n]*kmax*(Nk-1)/(Nt-1);
    p1 = exp(i*sum(k.*d1));
    p2 = exp(i*sum(k.*d2));
    p3 = exp(i*sum(k.*d3));
    p4 = exp(i*sum(k.*d4));
    g0 = (p1+p2+p3+p4)/4;
    g1 = (p1+p2-p3-p4)/4;
    g2 = (p1-p2+p3-p4)/4;
    g3 = (p1-p2-p3+p4)/4;
    h=[Esa/2 Vss*g0 0 0 0 Vsapc*g1 Vsapc*g2 Vsapc*g3 0 0;
       0 Esc/2 -Vpasc*conj(g1) -Vpasc*conj(g2) -Vpasc*conj(g3) 0 0 0 0 0;
       0 0 Epa/2 0 0 Vxx*g0 Vxy*g3 Vxy*g2 0 -Vpasec*g1;
       0 0 0 Epa/2 0 Vxy*g3 Vxx*g0 Vxy*g1 0 -Vpasec*g2;
       0 0 0 0 Epa/2 Vxy*g2 Vxy*g1 Vxx*g0 0 -Vpasec*g3;
       0 0 0 0 0 Epc/2 0 0 Vseapc*(g1) 0;
       0 0 0 0 0 0 Epc/2 0 Vseapc*(g2) 0;
       0 0 0 0 0 0 0 Epc/2 Vseapc*(g3) 0;
       0 0 0 0 0 0 0 0 Esea/2 0;
       0 0 0 0 0 0 0 0 0 Esec/2];
    %
    H = h + h';
    [V,D] = eig(H);
    eigst = sum(D);
    E(Nk,:) = sort(real(eigst));
    X(Nk) = -(Nk-1)/(Nt-1); %L-direction
    X1(Nk) = (Nk-1)/(Nt-1); %X-direction
end

hold on
h = plot(X,E,'b');
%h = plot(X1,E,'b');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
xlabel('k (as fraction of maximum value)')
ylabel('Energy (eV)')
grid on

%Note: X-axis from 0 to +1 represents the -X direction
%while the section from 0 to -1 represents the -L direction

