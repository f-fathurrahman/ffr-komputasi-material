import MyACE1x

r0 = 2.88
basis = MyACE1x.ace_basis(
    elements = [:Ti, :Al],
    order = 3,
    totaldegree = 6,
    rcut = 5.5,
    r0 = r0
);
