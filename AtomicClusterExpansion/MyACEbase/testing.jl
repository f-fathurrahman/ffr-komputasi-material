
module Testing

using Test, Printf

using ACEbase.FIO: read_dict, write_dict, save_dict, load_dict
using LinearAlgebra: norm
using StaticArrays

export print_tf, test_fio, h0, h1, h2, h3, println_slim 


function h0(str)
   dashes = "≡"^(length(str)+4)
   printstyled(dashes, color=:magenta); println()
   printstyled("  "*str*"  ", bold=true, color=:magenta); println()
   printstyled(dashes, color=:magenta); println()
end

function h1(str)
   dashes = "="^(length(str)+2)
   printstyled(dashes, color=:magenta); println()
   printstyled(" " * str * " ", bold=true, color=:magenta); println()
   printstyled(dashes, color=:magenta); println()
end

function h2(str)
   dashes = "-"^length(str)
   printstyled(dashes, color=:magenta); println()
   printstyled(str, bold=true, color=:magenta); println()
   printstyled(dashes, color=:magenta); println()
end

h3(str) = (printstyled(str, bold=true, color=:magenta); println())


print_tf(::Test.Pass) = printstyled("+", bold=true, color=:green)
print_tf(::Test.Fail) = printstyled("-", bold=true, color=:red)
print_tf(::Tuple{Test.Error,Bool}) = printstyled("x", bold=true, color=:magenta)

println_slim(::Test.Pass) = printstyled("Test Passed\n", bold=true, color=:green)
println_slim(::Test.Fail) = printstyled("Test Failed\n", bold=true, color=:red)

"""
`test_fio(obj): `  performs two tests:

- encodes `obj` as a Dict using `write_dict`, then decodes it using
`read_dict` and tests whether the two objects are equivalent using `==`
- writes `Dict` to file then reads it and decodes it and test the result is
again equivalent to `obj`

The two results are returned as Booleans.
"""
function test_fio(obj; warntype = true)
   D = write_dict(obj)
   test1 = (obj == read_dict(D))
   if !test1 
      @warn("test_fio (1) fails - read_dict(write_dict(obj)) != obj")
   end
   tmpf = tempname() * ".json"
   save_dict(tmpf, D)
   D2 = load_dict(tmpf)
   # if D != D2 
   #    @warn("load_dict(save_dict(D)) != D")
   # end
   obj2 = read_dict(load_dict(tmpf))
   if warntype && (typeof(obj) != typeof(obj2))
      @warn(
      """test_fio: the loaded object does not have the same type
             original : $(typeof(obj))
         deserialised : $(typeof(obj2))
      """)
   end
   test2 = (obj == obj2)
   if !test2 
      @warn("test_fio (2) fails - obj2 != obj")
   end
   return test1, test2
end

_Vec(X::AbstractVector{<: StaticVector{3}}) = 
      collect(reinterpret(Float64, X))

_svecs(x::AbstractVector{T}) where {T} = 
      collect(reinterpret(SVector{3, T}, x))


fdtest(F, dF, X::AbstractVector{<: StaticVector{3}}; kwargs...) = 
      fdtest( x -> F(_svecs(x)), 
              x -> _Vec(dF(_svecs(x))), 
              _Vec(X); kwargs... )

fdtest(F, dF, X::Number; kwargs...) = 
      fdtest( x -> F(x[1]), 
              x -> [dF(x[1])], 
              [X]; kwargs... )


"""
first-order finite-difference test for scalar F
```julia
fdtest(F, dF, x; h0 = 1.0, verbose=true)
```
"""
function fdtest(F, dF, x::AbstractVector; h0 = 1.0, verbose=true)
   errors = Float64[]
   E = F(x)
   dE = dF(x)
   # loop through finite-difference step-lengths
   verbose && @printf("---------|----------- \n")
   verbose && @printf("    h    | error \n")
   verbose && @printf("---------|----------- \n")
   for p = 2:11
      h = 0.1^p
      dEh = copy(dE)
      for n = 1:length(dE)
         x[n] += h
         dEh[n] = (F(x) - E) / h
         x[n] -= h
      end
      push!(errors, norm(dE - dEh, Inf))
      verbose && @printf(" %1.1e | %4.2e  \n", h, errors[end])
   end
   verbose && @printf("---------|----------- \n")
   if minimum(errors) <= 1e-3 * maximum(errors)
      verbose && println("passed")
      return true
   else
      @warn("""It seems the finite-difference test has failed, which indicates
      that there is an inconsistency between the function and gradient
      evaluation. Please double-check this manually / visually. (It is
      also possible that the function being tested is poorly scaled.)""")
      return false
   end
end

dirfdtest(F, dF, x, u; kwargs...) =
      fdtest(t -> F(x + t * u),
             t -> dF(x + t * u) .* Ref(u),
             0.0; kwargs...)


end
