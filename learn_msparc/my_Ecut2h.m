function h = my_Ecut2h(Ecut, FDn)
% Find mesh size h (Bohr) in finite-difference corresponding to a Ecut (Ha)
% in plane-wave codes.
% FDn: half of finite difference order

% this can be tuned to make the correspondence between h and Ecut more
% accurate
epsilon = 1e-1; % tolerance threshold for FD 2nd derivative approx.

% finite difference weights
w2 = zeros(1,FDn+1); 
for k=1:FDn
    w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
                    (k*k*factorial(FDn-k)*factorial(FDn+k));
    w2(1) = w2(1)-2*(1/(k*k));
end

kk = linspace(0,pi,1001);
y_cos =  -w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kk);
freq_err = abs(y_cos - kk.^2);
kc = kk(find(freq_err < epsilon,1,'last'));

h = kc / sqrt(2*Ecut);

end