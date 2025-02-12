% 0.5*Tr(X'*L*X) + alpha/4*rho(X)'*L^{dag}*rho(X)
function [f,g] = energy_func_01(X,~)
    LX = L*X;
    rhoX = sum(X.^2, 2); % diag(X*X');
    tempa = Lu\(Ll\rhoX); tempa = alpha*tempa;
    
    % function
    f = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
    % gradient
    g = LX + bsxfun(@times,tempa,X);
end