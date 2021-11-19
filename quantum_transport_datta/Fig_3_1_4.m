clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
epsil = 8.854e-12;
q = 1.602e-19;

%Lattice
Np = 200;
a = (10e-10/Np);
R = a*[1:1:Np];
t0 = (hbar^2)/(2*m*(a^2))/q;

% Hamiltonian,H = Kinetic,T + Potential,U + Uscf
T = (2*t0*diag(ones(1,Np))) - (t0*diag(ones(1,Np-1),1)) - (t0*diag(ones(1,Np-1),-1));
UN = (-q*2/(4*pi*epsil))./R; % Z=2 for Helium
Uscf = zeros(1,Np);
change = 1;
while change > 0.01
    [V,D] = eig(T + diag(UN+Uscf));
    D = diag(D);
    [DD,ind] = sort(D);
    E = D(ind(1));
    psi = V(:,ind(1));
    P = psi.*conj(psi);
    P = P';
    Unew = (q/(4*pi*epsil))*((sum(P./R) - cumsum(P./R)) + (cumsum(P)./R));
    change = sum(abs(Unew-Uscf))/Np;
    Uscf = Unew;
end

%analytical solutions for 1s hydrogen
a0 = 4*pi*epsil*hbar*hbar/(m*q*q);
P0 = (4*a/(a0^3))*R.*R.*exp(-2*R./a0);

hold on
%h=plot(R,UN, ' b ' );% Part (a)
%h=plot(R,Uscf, ' b ' );% Part(a)
h = plot(R, P, 'b');% Part (b)
h = plot(R, P0, 'bx');% Part (b)
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
xlabel('R (m)');
%ylabel('U (eV)'); % Part (a)
%axis([0 1e-9 -100 20]); % Part (a)
ylabel('Probability'); % Part (b)
axis([0 1e-9 0 0.1]); % Part (b)
grid on
