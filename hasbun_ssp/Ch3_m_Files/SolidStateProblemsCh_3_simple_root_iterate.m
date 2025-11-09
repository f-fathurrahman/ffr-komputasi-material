x=2.3;        %the initial guess
c=0.0120;     %the constant
for i=1:10    %the loop - 10 iterations do fine
x=-log(c/x^2) %the new guess, which will converge
end
