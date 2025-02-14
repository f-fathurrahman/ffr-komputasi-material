function g = PROB01_grad(X)
    rhoX = sum(X.^2, 2); % diag(X*X');
    tempa = Lu\(Ll\rhoX); tempa = alpha*tempa;
    g = L*X + bsxfun(@times,tempa,X);
end
