function [f, g] = PROB02_fun(X, Ham)
    XHX = X' * (Ham * X);
    [~, V] = eig(XHX);
    [Vs, ~] = sort(diag(V));
    f = trace(XHX);
    g = 2*( Ham*X - X*diag(Vs) );
end
