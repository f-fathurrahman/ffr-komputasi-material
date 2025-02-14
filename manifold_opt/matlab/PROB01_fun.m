function [f,g] = PROB01_fun(X, L, Lu, Ll, alpha)
    LX = L*X;
    rhoX = sum(X.^2, 2); % diag(X*X');
    tempa = Lu\(Ll\rhoX);
    tempa = alpha*tempa;
    f = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
    g = LX + bsxfun(@times,tempa,X);
end