clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
epsil = 8.854e-12;
q = 1.602e-19;

% Lattice
Np=200;
a=(10e-10/Np);
R=a*[1:1:Np];t0=(hbarˆ2)/(2*m*(aˆ2))/q;

% Hamiltonian,H = Kinetic,T + Potential,U + Ul + Uscf
T = (2*t0*diag(ones(1,Np)))-(t0*diag(ones(1,Np-1),1))-(t0*diag(ones(1,Np-1),-1));
UN = (-q*14/(4*pi*epsil))./R; % Z=14 for silicon
l = 1;
Ul = (l*(l+1)*hbar*hbar/(2*m*q))./(R.*R);
Uscf = zeros(1,Np);
change = 1;
while change>0.1
    [V,D] = eig(T + diag(UN+Uscf));
    D = diag(D);
    [DD, ind] = sort(D);
    E1s = D(ind(1));
    psi = V(:,ind(1));
    P1s = psi.*conj(psi);
    P1s = P1s';
    E2s = D(ind(2));
    psi = V(:,ind(2));
    P2s = psi.*conj(psi);
    P2s = P2s';
    E3s = D(ind(3));
    psi = V(:,ind(3));
    P3s = psi.*conj(psi);
    P3s = P3s';
    [V,D] = eig(T + diag(UN+Ul+Uscf));
    D = diag(D);
    [DD,ind] = sort(D);
    E2p = D(ind(1));
    psi = V(:,ind(1));
    P2p = psi.*conj(psi);P2p=P2p';
    E3p = D(ind(2));
    psi = V(:,ind(2));
    P3p = psi.*conj(psi);
    P3p = P3p';
    n0 = (2*(P1s + P2s + P3s)) + (6*P2p) + (2*P3p);
    n = n0*(13/14);
    Unew = (q/(4*pi*epsil))*((sum(n./R)-cumsum(n./R))+(cumsum(n)./R));
    %Uex=(-q/(4*pi*epsil))*((n./(4*pi*a*R.*R)).ˆ(1/3));%Unew=Unew+Uex;
    change = sum(abs(Unew-Uscf))/Np;
    Uscf = Unew;
end

[E1s E2s E2p E3s E3p]

% analytical solution for 1s hydrogen
a0 = 4*pi*epsil*hbar*hbar/(m*q*q);
P0 = (4*a/(a0ˆ3))*R.*R.*exp(-2*R./a0);
hold on
h = plot(R, P1s, 'b');
h = plot(R, P0, 'bx');
h = plot(R, P3p, 'bo');
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize' ,[25])
xlabel('R (m)');
ylabel('Probability');
axis([0 5e-10 0 0.08]);
grid on
