

#A = [1.0 0 2; 0 3 0; 4 0 5];
#A = sparse(A)
#x = [1.0, 2.0, 3.0]
#Ax = similar(x)

#@b my_CSC_mul!(Ax, A, x)

A = [1.0  0  0  2  0;
     0  3  0  0  0;
     4  0  5  0  6;
     0  7  0  8  0;
     9  0 10  0 11]
A = sparse(A)
x = collect(1:5) .+ 0.1
Ax = similar(x)

N = 10000;
A = sprand(ComplexF64, N, N, 0.01);
x = rand(ComplexF64, N);
Ax = similar(x);


