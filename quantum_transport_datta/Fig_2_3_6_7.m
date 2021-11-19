clear all
close all

% Constants (all MKS, except energy which is in eV)
hbar = 1.055e-34;
m = 9.110e-31;
epsil = 8.854e-12;
q = 1.602e-19;

a0 = 4*pi*epsil*hbar*hbar/(m*q*q);
E0 = q/(8*pi*epsil*a0);

% Lattice
Np =100;
a = (5e-10*2/Np); % *1 for Fig.1.3.6 and *2 for Fig.1.3.7
R = a*[1:1:Np];
t0 = (hbar^2)/(2*m*(a^2))/q;

%Quantum numbers
n = 1;
l = 0; % for 1s, n=1 and for 2s, n=2

%Hamiltonian, H = Kinetic,K + Potential,U
K = (2*t0*diag(ones(1,Np))) - (t0*diag(ones(1,Np-1),1)) - (t0*diag(ones(1,Np-1),-1));
U = ((-q/(4*pi*epsil)./R) + (l*(l+1)*hbar*hbar/(2*m*q))./(R.*R));
U = diag(U);

[V,D] = eig(K+U);
D = diag(D);
[DD, ind] = sort(D);
E = D(ind(n-l));
psi = V(:,ind(n-l));
P = psi.*conj(psi);

% Compare
[-E0/(n^2) E]

% analytical solutions
P1s = (4*a/(a0^3))*R.*R.*exp(-2*R./a0);
P2s = (4*a/(2*4*4*(a0^3)))*R.*R.*((2-(R./a0)).^2).*exp(-2*R./(2*a0));
P3s = (4*a/(3*81*81*(a0^3)))*R.*R.*((27-(18*R./a0)+(2*(R./a0).^2)).^2).*exp(-2*R./(3*a0));
P2p = (4*a/(3*32*(a0^3)))*R.*R.*((R./a0).^2).*exp(-2*R./(2*a0));
P3p = (8*a/(3*81*81*(a0^3)))*R.*R.*((6-(R./a0)).^2).*((R./a0).^2).*exp(-2*R./(3*a0));

hold on
h = plot(R, P, 'b');
h = plot(R, P1s, 'bx'); % use P1s for '1s' and P2s for '2s'
set(h, 'linewidth', [2.0])
set(gca, 'Fontsize', [25])
xlabel('x (m)');
ylabel('Probability');
grid on
