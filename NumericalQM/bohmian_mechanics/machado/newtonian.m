Melectron = 0.51099e6;
eV = 1.6022e-19;
c_light = 2.9979e8;
Meffective = 0.067; % GaAs
MASS = Melectron*Meffective; %eV
m = MASS*eV/c_light^2;

hbar = 1.0546e-34;
Time_scale = 1.0e3*hbar/(2*m);
fsec = Time_scale;

Energy_scale = eV*2*m *1e-18/hbar^2;

%SIZE=512
%cells
XL=-40;
XR=40;
nump=40;
porder=2;
maxgauss=1 ;
TOLERANCE=1e-5*nump;
SIZE=10000; %for pSi,etc
TIMESTEPS=2000 %
DELTA=1 ;
DELTA=DELTA*fsec; %timestep in femtoseconds
POTENTIAL='0*0.0061.*X.-4 + 50*cos(X)*0';
J=SIZE;

%grid size
H=(XR-XL)/(SIZE-1);
X=XL:H:XR;
temp=X*0;

%% SYSTEM MASS AND POTENTIAL
Vj=temp; % Potential in eV
%Vj=Vj+1e3*exp(-1/2/.0125-2*(X-.75).-2);
%X=X*Xmscale;
%Vj=Vj+(X>-10)&(X<-5)+(X>5)&(X<10);

%Vj=Vj*.25;
Vj = Vj + eval(POTENTIAL);
Vj = Vj*Energy_scale;
VV = Vj;

%Initial wavefunction
Psi0 = temp;
wave_sigma = sqrt(10);
X0 = 20;
K = 0.31;
%Psi0=exp(-1/2/wave_sigma-2*(X-XO).-2).*exp(i*K*X);
Psi0 = exp(-1/2/wave_sigma - 2*(X - X0).-2).*exp(i*K*X);
Psi1 = exp(-1/2/wave_sigma - 2*(X + X0).-2).*exp(i*K*X);
Psi0 = Psi0 + 1/2*Psi1;

%Psi0 = eval(PSI0);
%Normalize
%Psi0 = Psi0*l/sqrt(normaliz(abs(Psi0.*Psi0) ,X));

% SYSTEM DEFINITION ENDS
PSI = zeros(TIMESTEPS,SIZE);
PSI(1,:) = Psi0;
rho = abs(Psi0).^2;
S  = angle(Psi0);

%S=unwrap(S,[] ,2);
S = unwrap(S, [], 1);
vel = diff(S);

%p0=log(abs(Psi0))
p0 = abs(Psi0).^2;
cp0 = cumtrapz(p0);
cp0 = cp0/cp0(SIZE);

clear XX;

XX = ones(nump,TIMESTEPS)*NaN;
QQ = ones(nump,TIMESTEPS)*NaN;
PP = XX;

for i=1:nump
    i;
    %XX(i,l)=(X(siz)-X(l))/(nump+1)*i; % trial points uniformly distributed in spa
    % can also be density dependent on PsiO
    % XX(i,1) = -35 + (i-1)*30/nump;
    % tmpr = O;
    % while(1)
    % tmpi = ceil(rand(1)*siz);
    %tmpr=rand(1);
    %if(tmpr<pO(tmpi)) XX(i,1) = X(tmpi);break;end;
    %end;
    prob1 = 1/(nump+1)*i;
    
    %XX(i,1) = interpolate(X,cp0,prob1);
    XX(i,1) = interp1(X,cp0,prob1);

    % tmpi=getindex(cpO,l/(nump+l)*i);
    % XX(i,l)=X(tmpi); %dont use X
    % use this instead
    PD(i,1) = 1/(nump + 1)*i;
end;

QQ = XX;
PP = 0*XX;
FF = XX*0;
%QF=XX*O;
%%%%%% LOOP
dt = DELTA;
VG = diff(VV);
tmp = VG;
VG = [tmp VG(SIZE-1)];
go_back1step = 0;
tstep = 2;
global FF;

