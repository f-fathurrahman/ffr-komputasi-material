function npl = my_mesh2chebdegree(h) 
    % the relation between h and npl is fit with a cubic polynomial
    % p(x) = p3 * x^3 + p2 * x^2 + p1 * x + p0.
    p3 = -700. / 3.;
    p2 = 1240. / 3.;
    p1 = -773. / 3.;
    p0 = 1078. / 15.;
    if (h > 0.7) 
        npl = 14;
    else 
        npl = ((p3 * h + p2) * h + p1) * h + p0;
    end
    npl = round(npl);
end