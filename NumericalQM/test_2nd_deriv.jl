import PyPlot
const plt = PyPlot

include("init_FD1d_grid.jl")
include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")

function my_gaussian(x, α=1.0)
    return exp(-α*x^2)
end

function d2_my_gaussian(x, α=1.0)
    return -2*α*exp(-α*x^2) + 4*α^2 * x^2 * exp(-α*x^2)
end

function main()
    N = 50
    x, h = init_FD1d_grid( (-5.0, 5.0), N )
    fx = my_gaussian.(x)

    plt.clf()
    plt.plot(x, fx, marker="o")
    plt.grid()
    plt.savefig("TEMP_my_gaussian.pdf")

    d2_fx = d2_my_gaussian.(x)

    D2_3pt = build_D2_matrix_3pt(N, h)
    D2_fx_3pt = D2_3pt*fx

    D2_5pt = build_D2_matrix_5pt(N, h)
    D2_fx_5pt = D2_5pt*fx

    plt.clf()
    plt.plot(x, d2_fx, label="analytic")
    plt.plot(x, D2_fx_3pt, label="FD-3pt")    
    plt.plot(x, D2_fx_5pt, label="FD-5pt")
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_d2_my_gaussian.pdf")


    plt.clf()
    plt.plot(x, D2_fx_3pt - d2_fx, label="FD-3pt")    
    plt.plot(x, D2_fx_5pt - d2_fx, label="FD-5pt")
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_d2_my_gaussian_diff.pdf")
end

main()
