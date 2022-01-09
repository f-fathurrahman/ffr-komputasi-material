clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
q = 1.602e-19;
a = 3e-10;
m = 9.110e-31;

% Conduction band parameters
mw = 0.07*m;
ma = 0.22*m;
mb = 0.7*mw + 0.3*ma;
kk = 0*.1*pi;
Ec = 0;
Eb = (.7*0) + (.3*1.25);
for nk = 1:24
    Nw = nk + 10;
    Nb = 2*Nw;
    Np = Nb + Nw + Nb;
    W(nk) = (Nw-1)*a*1e9;
    tb = (hbar^2)/(2*mb*(a^2)*q);
    tw = (hbar^2)/(2*mw*(a^2)*q);
    t = [tb*ones(1,Nb) tw*ones(1,Nw-1) tb*ones(1,Nb)];
    tt = [0 t] + [t 0];
    Ebk = Eb + (tb*(kk^2));
    Ewk = tw*(kk^2);
    Ebwk = (Eb/2)+((tb+tw)*(kk^2)/2);
    U = Ec + [Ebk*ones(1,Nb) Ebwk Ewk*ones(1,Nw-2) Ebwk Ebk*ones(1,Nb)];
    H = -diag(t,1) - diag(t,-1) + diag(tt) + diag(U);
    [V,D] = eig(H);
    D = diag(D);
    D = (sort(real(D)))';
    E1(nk) = D(1);
    E2(nk) = D(2);
end

hold on
h1 = plot(W, E1,' b');
h1 = plot(W, E2,' b--');
set(h1, 'linewidth',[2.0])
set(gca, 'Fontsize',[24])
xlabel('W (nm)')
ylabel('Energy (eV)')
axis([2 10 0 .4])
grid on
