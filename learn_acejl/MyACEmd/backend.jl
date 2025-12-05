
@inline function ace_evaluate!(tmp, calc, R::AbstractVector, species::AbstractVector, species0 )
    return MyACE1.evaluate!(tmp, calc, R, species, species0)
end

@inline function ace_evaluate!(B, tmp, calc, R::AbstractVector, species::AbstractVector, species0 )
    MyACE1.evaluate!(B, tmp, calc, R, species, species0)
    return B
end

function ace_evaluate(calc, R::AbstractVector, species::AbstractVector, species0)
    tmp = MyACE1.alloc_temp(calc, length(R))
    return ace_evaluate!(tmp, calc, R, species, species0)
end

function ace_evaluate(calc::MyACE1.IPBasis, R::AbstractVector, species::AbstractVector, species0)
    B = MyACE1.alloc_B(calc)
    tmp = MyACE1.alloc_temp(calc, length(R))
    return ace_evaluate!(B, tmp, calc, R, species, species0)
end


@inline function ace_evaluate_d!(out, tmp, calc, R::AbstractVector, species::AbstractVector, species0)
    e = MyACE1.evaluate_d!(out, tmp, calc, R, species, species0 )
    return e, tmp
end

function ace_evaluate_d(calc, R::AbstractVector, species::AbstractVector, species0)
    tmp = MyACE1.alloc_temp_d(calc, length(R))
    return ace_evaluate_d!(tmp.dV, tmp, calc, R::AbstractVector, species::AbstractVector, species0)
end


function ace_evaluate_d(basis::MyACE1.IPBasis, R::AbstractVector, species::AbstractVector, species0)
    dB  = MyACE1.alloc_dB(basis, length(R))
    tmp = MyACE1.alloc_temp_d(basis, length(R))
    ace_evaluate_d!(dB, tmp, basis, R, species, species0)
    return dB
end