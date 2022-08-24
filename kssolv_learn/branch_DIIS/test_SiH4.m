clear all;
close all;
randn('state', 0);
rand('state', 0);
cd ../;
KSSOLV_startup;
cd example;
sih4_setup;
%si8_setup;
%si64_setup;
cd ../src/GW;
ksinfo = gwsetup(mol);
omega = 0;
eta = 0;
[ng, nr] = size(ksinfo.F);
nclimit = 300;
chi0 = getchi0(ksinfo, omega, eta, nclimit);
epsilon = geteps(chi0, ksinfo.coulG);
inveps = inv(epsilon);
W = getw(ksinfo, chi0);
nvbands = 4;
ncbands = 1;
%[Esx, Ech, Ex] = getgw(ksinfo, eta, nvbands, ncbands, chi0);
%diagEsx = diag(Esx)*27.2116
%diagEch = diag(Ech)*27.2116
%diagEx = diag(Ex)*27.2116
%diagEx = diag(Ex)
[d, VA, WA, VB, WB, Hbse] = getkernel(ksinfo, omega, eta, nvbands, ncbands, W);
d
Hbse
[X1, X2, Lambda] = BSE_complex(diag(40.0*d)+Hbse, VB+WB);
lambda1 = diag(Lambda)
save('../../lambda1.mat','lambda1');
vcrank_ratio = 1.00;
vvrank_ratio = 1.00;
ccrank_ratio = 1.00;
[d, VA, WA, VB, WB, Hbse] = getkernelisdf(ksinfo, omega, eta, nvbands, ncbands, W, vcrank_ratio, vvrank_ratio, ccrank_ratio);
d
Hbse
[X1, X2, Lambda] = BSE_complex(diag(40.0*d)+Hbse, VB+WB);
lambda2 = diag(Lambda)
errorLambda = norm(lambda1 - lambda2)/norm(lambda1)
errorMax = max(lambda1 - lambda2)
cd ../../;
