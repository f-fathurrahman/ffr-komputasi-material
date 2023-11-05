module Eigensolvers

using Test
using LinearAlgebra
using Plots

# A,B and P must implement *
# A must be hermitian, A implemented as a function (not a matrix)
# B and P must be HDP
Id(x) = x
function eig_pcg(A,x0;B=Id,P=Id,tol=1e-16,maxiter=400,do_cg=true)
    #pre-allocate big arrays
    X = zeros(eltype(x0),size(x0)[1],3)
    AX = zeros(eltype(x0),size(x0)[1],3)
    BX = zeros(eltype(x0),size(x0)[1],3)

    # ensure normalization of x0
    Bx0 = B(x0)
    x = (1 / sqrt(real(dot(x0,Bx0))))*x0
    Ax = A(x)
    Bx = Bx0 / sqrt(real(dot(x0,Bx0)))
    @test isapprox(dot(x,Bx),1)
    lambda = 0
    x_prev = zeros(eltype(x0),size(x0)[1])
    Ax_prev = zeros(eltype(x0),size(x0)[1])
    Bx_prev = zeros(eltype(x0),size(x0)[1])
    res_history = zeros(maxiter)

    for i=1:maxiter
        lambda_prev = lambda
        # build new (preconditioned) residual and convergence test
        lambda = dot(x,Ax) / dot(x,Bx)
        @test abs(imag(lambda)) < 1e-9
        lambda = real(lambda)
        resid = Ax - lambda*Bx
        println("iter $i, λ = $lambda, res = $(norm(resid))")
        res_history[i] = norm(resid)
        if norm(resid) < sqrt(tol) || abs(lambda - lambda_prev) < tol
            return x,lambda, res_history[1:i]
        end
        # P_resid = P(resid)
        # AP_resid = A(P_resid)
        # BP_resid = B(P_resid)

        nsub = (do_cg && (i > 1)) ? 3 : 2

        # build subspace
        X[:,1] = x
        AX[:,1]= Ax
        BX[:,1]= Bx
        X[:,2] = P(resid)
        AX[:,2] = A(X[:,2])
        BX[:,2] = B(X[:,2])
        if nsub == 3
            X[:,3] .= x.-x_prev
            AX[:,3] .= Ax.-Ax_prev
            BX[:,3] .= Bx-Bx_prev
        end

        # bookkeeping previous x
        x_prev = x
        Ax_prev = Ax
        Bx_prev = Bx

        # rotate into new x
        @views begin
            c = rayleigh_ritz(X[:,1:nsub],AX[:,1:nsub],BX[:,1:nsub])
            x = X[:,1:nsub]*c
            Ax = AX[:,1:nsub]*c
            Bx = BX[:,1:nsub]*c
        end
        @test isapprox(dot(x,Bx),1)
    end

    return x,lambda, res_history
end

# Returns a column vector c such that X c is B-normalized and achieves the minimum of the Rayleigh quotient in the subspace spanned by X
function rayleigh_ritz(X, AX, BX)
    hamiltonian = X'*AX
    overlap = X'*BX
    @test isapprox(hamiltonian,hamiltonian',atol=1e-9)
    @test isapprox(overlap, overlap',atol=1e-9)
    hamiltonian = (hamiltonian + hamiltonian')/2
    overlap = (overlap + overlap')/2

    # ## Simpler and faster code, but fails if overlap is badly conditioned
    # D,V = eig(hamiltonian,overlap)
    # c = V[:,1]


    #TOFIX this is still not good enough, and accuracy saturates at 1e-7...
    #A x = lambda B x. Assume B = V S V^* = (V S^1/2) (V S^1/2)^*, then the problem is equivalent to S^-1/2 V^* A V S^-1/2 y = lambda y, with y = S^1/2 V^* x. Here we also allow dropping columns of B in case of ill-conditioning
    F = eigen(overlap)
    n = count(x -> (x < 1e-14),F.values)
    V = F.vectors[:,n+1:end]
    S = diagm(0=>F.values[n+1:end]) # could be optimized
    ham_reduced = S^(-1/2)*V'*hamiltonian*V*S^(-1/2)
    ham_reduced = (ham_reduced + ham_reduced')/2
    c = eigvecs(ham_reduced)[:,1]
    @test isapprox(norm(c),1)
    c = V*(S^(-1/2)*c)

    @test isapprox(dot(c,overlap*c),1)
    c = c / sqrt(real(dot(c,overlap*c))) # ensure normalization, normally not needed if eigensolver recognizes ham and ovl are hermitian
    return c
end

function test_pcg()
    seed!(0)
    n = 10
    A = randn(n,n) + im*randn(n,n)
    A = (A+A')/2 # force hermitian
    B = randn(n,n) + im*randn(n,n)
    B = (B'*B) + I # force HDP
    x0 = randn(n)
    x,lambda, res_history = eig_pcg(x->A*x,x0,B=(x->B*x))
    println(lambda,eigvals(A,B)[1])
    @test isapprox(lambda,eigvals(A,B)[1],rtol=1e-6)
    plot(res_history, yaxis=:log10)
end


"""
conjugate-gradient algorithm, works for positive definite matrix
"""
function cg(A,b;P=Id,x0=b,tol=1e-10,Imax=1000)
   #Initialization
   i=0
   g = A(x0)-b
   h = -P(g)
   Pg = -h
   x=x0 #need to pre-allocate x otherwise x is only local in the while-loop

   #iteration
   while i < Imax && norm(h) > tol
      i +=1
      Ah = A(h) #store matrix-vector product to minimize operations
      rho = - g'*h/(h'*Ah)
      x += rho*h
      g_prev = g
      g += rho*Ah
      Pg_prev = Pg #store
      Pg = P(g)
      gam = g'*Pg/(g_prev'*Pg_prev)
      h = -Pg + gam*h
   end
   return x, i, norm(h)
end

"""
eig_lanczos returns the eigenvalues of Ax = λBx where B is positive-definite
use selective orthogonalisation (do_so=true) if numerical instabilities 
"""

function eig_lanczos(A,x0; m=3, B=Id, tol=1e-5, Imax=min(1000,size(x0)[1]), norm_A = norm(A(x0))/norm(x0), do_so = false)
   #output
   vp = zeros(m)
   Vp = zeros(eltype(x0),size(x0)[1],m)
   resid = zeros(Imax) #taille des résidus
   resid[1] = 1.

   #Lanczos matrices
   diag_T = zeros(Imax)
   off_diag_T = zeros(Imax-1) #T = SymTridiagonal(diag_T,off_diag_T)
   V = zeros(eltype(x0),size(x0)[1],Imax) #contains the bi-orthogonal vectors v_j
   W = zeros(eltype(x0),size(x0)[1],Imax) #contains the bi-orthogonal vectors w_j

   Q = Matrix{eltype(x0)}(I,Imax,Imax) #stores v_j* B v_k for selective orthogonalization

   if do_so
      for k in 1:Imax-1
         Q[k,k+1] = eps()
         Q[k+1,k] = eps()
      end
   end

   #Initialization
   i=1
   w_temp = B(x0)
   b = sqrt(x0'*w_temp)
   @test abs(real(b)) > 1e-14
   W[:,1] = 1/b*w_temp
   V[:,1] = 1/b*x0
   r_temp = A(V[:,1])
   a=V[:,1]'*r_temp
   @test abs(imag(a)) < 1e-10
   diag_T[1] = real(a)
   r = r_temp - diag_T[1]*W[:,1]
   q = cg(B,r,x0=r)[1] #tol=1e-10, Imax=1000

   #iteration
   while resid[i] > tol && i < Imax
      i += 1
      b = sqrt(q'*r)
      @test abs(imag(b)) < 1e-10
      off_diag_T[i-1] = real(b)
      @test off_diag_T[i-1] > 1e-14
      V[:,i]=q/off_diag_T[i-1]
      W[:,i] = r/off_diag_T[i-1]

      if do_so && i>2
         omega_temp = off_diag_T[1]*Q[i-1,2] + (diag_T[1]-diag_T[i-1])*Q[i-1,1] - off_diag_T[i-2]*Q[i-2,1]
         Q[i,1] = (omega_temp + 2sign(omega_temp)*norm_A*eps())/off_diag_T[1]
         for k in 2:i-1
            omega_temp = off_diag_T[k]*Q[i-1,k+1] + (diag_T[k]-diag_T[i-1])*Q[i-1,k] + off_diag_T[k-1]*Q[i-1,k-1] - off_diag_T[i-2]*Q[i-2,k]
            Q[i,k] = (omega_temp + 2sign(omega_temp)*norm_A*eps())/off_diag_T[i-1]
         end
         #re-orthogonalisation condition
         if norm(Q[i,1:i-1],Inf) > sqrt(eps())
            println("Iteration : $(i), Off-set : $(norm(Q[i,1:i-1],Inf))")
            for k in 1:i-2
               Q[i-1,k] = eps() #re-orthogonalize => reset values to eps()
               Q[i,k] = eps()
               V[:,i-1] -= (W[:,k]'*V[:,i-1])*V[:,k]
               V[:,i] -= (W[:,k]'*V[:,i])*V[:,k]
            end
            W[:,i-1] = B(V[:,i-1])
            Q[i,i-1] = eps()
            V[:,i] -= (W[:,i-1]'*V[:,i])*V[:,i-1]
            W[:,i] = B(V[:,i])
         end
      end

      r_temp = A(V[:,i])
      a = r_temp'*V[:,i]
      @test abs(imag(a)) < 1e-10
      diag_T[i] = real(a)
      r = r_temp - diag_T[i]*W[:,i] - off_diag_T[i-1]*W[:,i-1]
      vp, Vp = eigen(SymTridiagonal(diag_T[1:i],off_diag_T[1:i-1]))
      for k in 1:min(i,m)
         resid[i] += off_diag_T[i-1]*abs(Vp[i,k])/m
      end
      println("iter $i, λ = $(vp[1]), res = $(resid[i])")
      q,i_cg,cg_tol = cg(B,r,x0=r) #tol=1e-10, Imax=1000
      @test cg_tol < 1e-8
      #q = B\r
   end
   return (V[:,1:i]*Vp)[:,1:min(i,m)], vp[1:min(i,m)], resid
end

end # module
