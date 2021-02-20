using Printf

include("root_brent.jl")

function f(x)
    return sin(sqrt(x)) - x
end

function main()
    xroot = root_brent(f, 0.5, 1.0)
end

main()
