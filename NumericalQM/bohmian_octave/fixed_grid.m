%% SYSTEM TIME AND SPACE
SIZE = 512
%cells
XL = -40.0
XR = 40.0; %nanometers

TIMESTEPS = 2000

eV   = 1.6022e-19;
hbar = 1.0546e-34;
c_light = 2.9979e8;
Melectron = 0.51099e6;

Meffective = 0.067; % GaAs
MASS  = Melectron*Meffective; %eV
MASSE = MASS;
MASS_KG = MASS*eV /c_light^2;

m = MASS*eV/c_light^2;
MASS = 1/2;
Lenght_scale = 1;
Time_scale = 1e3*hbar/(2*m);
Energy_scale = eV*2*m *1e-18/hbar^2;
Velocity_scale = 1e-3 * 2 * m/hbar;
M = m;
fsec = Time_scale;

%Check
cprime = c_light/1e6*Velocity_scale;

DELTA = 0.1 * fsec;

%timestep in femtoseconds
J = SIZE;
H = (XR-XL)/(SIZE-1);

%grid size
X    = XL:H:XR;
temp = X*0;

%% SYSTEM MASS AND POTENTIAL
Vj = temp; % Potential in eV

% Vj=Vj+ 1e3*exp{ -1/2/ .0125" 2*(X-. 75). "2);
%X=X*Xmsca Ie;
Vj = Vj + (X>-10) & (X<-5) + (X>5) & (X<10);
Vj = Vj*0.25;
Vj = Vj*Energy_scale;


%Initial wavefunction
Psi0 = temp;
wave_sigma = sqrt(10.0);
X0 = -25;
K = 0.31;
Psi0 = exp( -1/2/wave_sigma^2*(X-X0).^2) .* exp(i*K*X);

% Normalize
% PsiO=PsiO*1 /sqrt( normaliz( abs(PsiO. *PsiO),X));
% SYSTEM DEFINITION ENDS


PSI = zeros(TIMESTEPS,SIZE);
PSI(1,:) = Psi0;
ALPHAj = temp;
%Fixed
Gj=Vj-2*i/DELTA;
Gj=Gj*(2*MASS);

Dj = 1 - H^2/12*Gj;
Aj = 1 + H^2*Gj./(2*Dj);
j = J;
niuJ = (1 - abs( Aj(j)^2))/abs(1 - Aj(j)^2);

j = 1;
niu1 = (1 - abs( Aj(j)^2))/abs(1 - Aj(j)^2);

% If V1=VJ then niu1=niuJ symmetrical
% Calculate all LEGENDRE coeficients (#timesteps) ahead of time index
% 1->n=0
%REPEAT 1 For X=1
j = 1;
lambda = 2*H^2/DELTA; %Fixed
c = 1 - i*lambda/( 6*Dj(j) ); %Fixed

phi = angle((Aj(j).^2 - 1) ./ c); %Fixed

niu = (1 - abs(Aj(j)^2)) / abs(1-Aj(j)^2);

% PLEASE NOTE PN(n) is actually Legrendre polynomial of order n-1
% So use to refer to PN(m) m=n+1
PN = temp;
LN = temp;
PN(1) = 1;
PN(2) = niu; % These are PNO and PN1

for m = 2:(TIMESTEPS-1)
    %generate Legendre polynomials;
    n = m-1;
    PN(m+1) = 2*niu*PN(m) - PN(m-1) - (niu*PN(m) - PN(m-1))/(n+1);
end;

LN(1) = niu*exp(-i*phi); % LN(O) would be -1; LNn is LN(n)
for n = 2:(TIMESTEPS-1)
    m = n+1;
    LN(n) = exp(-i*n*phi)/(2*n-1)*(PN(m) - PN(m-2));
    %LN(n)=exp(-i*n*phi)*( niu*PN(m-1 )-PN( m-2) )*(2*n+ 1)/ (2*n-1)/ (n+ 1);
end;
LN1 = LN;

%REPEAT 2 for X=J
j=J;

lambda = 2*H^2/DELTA; %Fixed

c = 1 - i*lambda/(6*Dj(j)); %Fixed

phi = angle((Aj(j).^2-1)./c); %Fixed

niu = (1 - abs(Aj(j)^2))/abs(1-Aj(j)^2);

% PLEASE NOTE PN(n) is actually Lagrage polynomial of order n-1
% So use to refer to PN(m) m=n+1
PN = temp;
LN = temp;
PN(1) = 1;
PN(2) = niu; % These are PNO and PN1
for m = 2:(TIMESTEPS-1)
    %generate Legendre polynomials;
    n = m-1;
    PN(m+1) = 2*niu*PN(m) - PN(m-1) - (niu*PN(m) - PN(m-1))/(n + 1);
end;
LN(1) = niu*exp(-i*phi); % LN(O) would be -1; LNn is LN(n)

for n=2:(TIMESTEPS-1)
    m = n+1;
    LN(n)=exp(-i*n*phi) /(2*n-1)*(PN(m)-PN(m-2));
    %LN( n )=exp( -i*n*phi)*( niu*PN( m-1 )-PN( m-2) )*(2*n+1) /(2*n-1) /( n+ 1);
end;
LNJ = LN;

%###################################
%###################################
% LOOP
%###################################
%###################################
for N_INDEX=1:(TIMESTEPS-1)

    N_INDEX
    %###################################
    %STEP 1
    %###################################
    Fj = 4*i*PSI(N_INDEX,:)/DELTA;
    Fj = Fj*2*MASS;
    
    %###################################
    %STEP 2
    %###################################
    % From the left
    % Calc Ej Fixed
    t1 = Aj(1) + sqrt(Aj(1)^2-1);
    t2 = Aj(1) - sqrt(Aj(1)^2-1);
    if(abs(t1) > 1)
        ALPHAj(1) = t1;
    else
        ALPHAj(1) = t2;
    end;
    Ej(1) = ALPHAj(1);
    %Recurrence
    for j=2:SIZE
        Ej(j) = 2*Aj(j) - 1.0/Ej(j-1);
    end

    % Calc Qj Variable
    j = 1;
    %convolution
    SUM=0;
    for k=1:N_INDEX
        n = N_INDEX - k + 1;
        SUM = SUM + PSI(k,j)*LN1(n);
    end;
    j = 1;
    
    Qj(j) = conj(Dj(j)) * ( conj(Aj(j)) - ALPHAj(j) ) * PSI(N_INDEX,j) + Dj(j)*( Aj(j) - ALPHAj(j) )*SUM;

    %Recurrence
    for j=2:SIZE
        Qj(j) = Qj(j-1)/Ej(j-1) + H^2*Fj(j)/Dj(j);
    end
    
    %###################################
    %STEP 3
    %###################################
    %From the right
    t1 = Aj(J) + sqrt( Aj(J)^2 - 1 );
    t2 = Aj(J) - sqrt( Aj(J)^2 - 1);
    if (abs(t1) > 1)
        ALPHAj(J) = t1;
    else
        ALPHAj(J) = t2;
    end;
    j = J;
    %convolution
    SUM = 0;
    for k = 1:N_INDEX
        n = N_INDEX - k + 1;
        SUM = SUM + PSI(k,j)*LNJ(n);
    end;
    j = J;

    BETAj(j) = conj(Dj(j))*( conj(Aj(j)) - ALPHAj(j) )*PSI(N_INDEX,j) + Dj(j)*( Aj(j) - ALPHAj(j))*SUM;
    Wj(J) = ( Qj(J-1) + BETAj(J)*Ej(J-1) ) / ( 1 - ALPHAj(J)*Ej(J-1) );
    %Recurrence for Wj
    for j=(J):-1:2
        Wj(j-1) = ( Wj(j) - Qj(j-1)) /Ej(j-1);
    end

    %##############################
    %STEP 4
    %##############################
    % Calculate Psi from Wj
    PSI(N_INDEX+1,:) = PSI(N_INDEX,:) .* ( i*H^2./(3*DELTA*Dj) - 1 ) + Wj./Dj;
end;

figure
image( abs(PSI)*100 );
ylabel('Timestep (.1 fs)');
xlabel('X From -40 to 40 (nm)');


