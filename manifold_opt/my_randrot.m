% Source: https://github.com/sgmanohar/matlib
function Q = my_randrot(N)
% generate a random rotation matrix 
% Q = randrot(N) 
% N: dimensions 
% sgm 2019
Q = qr(randn(N)); 
