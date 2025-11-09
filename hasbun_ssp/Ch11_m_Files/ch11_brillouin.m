%copyright by J. E Hasbun and T. Datta
% ch11_brillouin.m

% This script plots the Brillouin function for different values of
% the total angular momentum J. The paramagnetic magnetization follows
% the Brillouin function curve.

close all; clear all;

% Constants
mub = 9.27*10^-24; kb = 1.38*10^-23;

% defining x = B/T (magnetic field/temperature ratio)

% Function definitions
landeg = @(L,S) 1.5 + (S*(S+1) - L*(L+1))/(2*(L+S)*(L+S+1));
brillouinJ = @(x,L,S)((2*(L+S)+1)/(2*(L+S)))*coth(((2*(L+S)+1)/...
    (2*(L+S)))*landeg(L,S)*(mub/kb)*(L+S)*x)-(1/(2*(L+S)))*coth((landeg(L,S)*...
    (mub/kb)*(L+S)*x)/(2*(L+S)));

% Printing Landeg values
fprintf('landeg(0.0,0.5)= %0.3f\n',landeg(0.0,0.5));
fprintf('landeg(1.0,0.5)= %0.3f\n',landeg(1.0,0.5));
fprintf('landeg(2.0,0.5)= %0.3f\n',landeg(2.0,0.5));
fprintf('landeg(3.0,0.5)= %0.3f\n',landeg(3.0,0.5));
fprintf('landeg(3.0,3.5)= %0.3f\n',landeg(3.0,3.5));
fprintf('landeg(6.0,2.0)= %0.3f\n',landeg(6.0,2.0));
fprintf('landeg(5.0,10.0)= %0.3f\n',landeg(5.0,10.0));


% Plotting the Brillouin function curves for different J values
xmax = 4;
hold on;
fplot(@(x)brillouinJ(x,0.0,0.5),[0,4],'-k', 'LineWidth',2);
fplot(@(x)brillouinJ(x,1.0,0.5),[0,4],'--r');
fplot(@(x)brillouinJ(x,2.0,0.5),[0,4],'-g');
fplot(@(x)brillouinJ(x,3.0,0.5),[0,4],'--b');
fplot(@(x)brillouinJ(x,3.0,3.5),[0,4],'-c');
fplot(@(x)brillouinJ(x,6.0,2.0),[0,4],'--m');
fplot(@(x)brillouinJ(x,5.0,10.0),[0,4],'-k', 'LineWidth',5);
xlabel('B/T');ylabel('Brillouin Function');
legend('J=0.5','J=1.5','J=2.5','J=3.5','J=6.5','J=8.0','J=15.0');
hold off;
