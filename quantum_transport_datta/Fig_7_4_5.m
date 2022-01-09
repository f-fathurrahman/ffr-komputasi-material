clear all
close all

%Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
a = 3e-10;
m = 9.110e-31;

%Luttinger-Kohn parameters
g1 = 6.85; g2 = 2.1; g3 = 2.9;%GaAs
w1 = (hbar^2)*g1/(2*m*q*(a^2));
w2 = (hbar^2)*g2/(2*m*q*(a^2));
w3 = (hbar^2)*g3/(2*m*q*(a^2));

g1 = 3.45; g2=0.68; g3=1.29;%AlAs
a1 = (hbar^2)*g1/(2*m*q*(a^2)); b1=(0.7*w1) + (0.3*a1);
a2 = (hbar^2)*g2/(2*m*q*(a^2)); b2=(0.7*w2) + (0.3*a2);
a3 = (hbar^2)*g3/(2*m*q*(a^2)); b3=(0.7*w3) + (0.3*a3);

Ev = 0;
Evb = (0.7*0) + (0.3*0.75)

Nw=18;
Nb = Nw;
Np = Nb + Nw + Nb;
W = (Nw-1)*a*1e9
Z = zeros(Np);

for nk=1:26
    k(nk) = (nk-1)/500;% in A^-1
    l = 0;
    m = 1;
    lm = sqrt((l^2)+(m^2));
    kx = (l/lm)*k(nk)*a*1e10;
    ky = (m/lm)*k(nk)*a*1e10;
    k2 = (kx^2) + (ky^2);

    t = [b1*ones(1,Nb) w1*ones(1,Nw-1) b1*ones(1,Nb)];
    tt = [0 t] + [t 0];
    Ebk = Evb + (b1*k2);
    Ewk = (w1*k2);
    Ebwk = (Ebk + Ewk)/2;
    U = Ev + [Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    P = -diag(t,1) - diag(t,-1) + diag(tt) + diag(U);

    t = -2*[b2*ones(1,Nb) w2*ones(1,Nw-1) b2*ones(1,Nb)];tt=[0 t]+[t 0];
    Ebk = b2*k2;
    Ewk = w2*k2;
    Ebwk = (Ebk+Ewk)/2;
    U = [Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    Q = -diag(t,1)-diag(t,-1)+diag(tt)+diag(U);

    Ebk=-(sqrt(3)*b2*((kx^2)-(ky^2)))+(i*2*b3*sqrt(3)*kx*ky);
    Ewk=-(sqrt(3)*w2*((kx^2)-(ky^2)))+(i*2*w3*sqrt(3)*kx*ky);
    Ebwk=(Ebk+Ewk)/2;
    U=[Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    R=diag(U);

    t=-2*i*sqrt(3)*(kx-(i*ky))*[b3*ones(1,Nb) w3*ones(1,Nw-1) b3*ones(1,Nb)]/2;
    S=diag(t,1)-diag(t,-1);

    H = [P+Q Z; Z P+Q];
    HL = [P-Q Z; Z P-Q];
    HC = [-S R; R' S'];    
    H = -[H HC;HC' HL];
    %[nk sum(sum(abs(H-H')))]

    [V,D]=eig(H);D=diag(D);D=-(sort(real(-D)))';
    E1(nk)=D(1);E2(nk)=D(2);E3(nk)=D(3);E4(nk)=D(4);
end

k = k*10;%per Angstrom to per nm
hold on
%h=plot(W,Ean1,'b');
%h=plot(W,Ean2,'b');
h = plot(k, E1, 'b');
set(h, 'linewidth', [2.0])
%h=plot(k,E2,'bx');
h = plot(k, E3, 'b');
%h=plot(k,E4,'b+');
set(h,'linewidth', [2.0])
set(gca,'Fontsize', [24])
xlabel('k (/nm)')
ylabel('Energy (eV)')
axis([0 0.5 -0.1 0])
grid on
