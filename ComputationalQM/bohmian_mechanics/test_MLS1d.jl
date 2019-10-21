using Printf
using LinearAlgebra
import PyPlot

const plt = PyPlot

include("MLS1d_shape.jl")

function test_main()

    l = 10.0
    dx = 0.1

    # Nodal points
    xi = 0.0:dx:l
    Nnodes = length(xi)

    # Coordinates of evaluation
    x = 0.0:0.05:l
    Npoints = length(x)

    # DETERMINE RADIUS OF SUPPORT OF EVERY NODE
    scl = 3
    dm = scl * dx * ones(Nnodes)

    # Evaluate MLS shape function at all evaluation points x
    ϕ, dϕ, d2ϕ = MLS1d_shape(3, Nnodes, xi, Npoints, x, dm, "gauss", 3.0)

    println("Done doing MLS, now printing")
    plt.clf()
    for i in 1:Nnodes
        plt.plot(x, ϕ[:,i], marker="o")
    end
    plt.savefig("TEMP_shape_0.png", dpi=150)

    plt.clf()
    for i in 1:Nnodes
        plt.plot(x, dϕ[:,i], marker="o")
    end
    plt.savefig("TEMP_shape_1.png", dpi=150)

    plt.clf()
    for i in 1:Nnodes
        plt.plot(x, d2ϕ[:,i], marker="o")
    end
    plt.savefig("TEMP_shape_2.png", dpi=150)


    yi  = sin.(xi)   # Nodal function values
    y   = sin.(x)   # Exact solution
    yh  = ϕ * yi  # Approximation function
    err = norm(y - yh) / norm(y) * 100  # Relative error norm in approximation function

    println("err (percent) = ", err)

    plt.clf()
    plt.plot(x, y, label="exact")
    plt.plot(x, yh, label="approx")
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_approx_sin.png", dpi=150)


    dy   = cos.(x)
    dyh  = dϕ * yi
    err = norm(dy - dyh) / norm(dy) * 100

    println("diff err (percent) = ", err)

    plt.clf()
    plt.plot(x, dy, label="exact")
    plt.plot(x, dyh, label="approx")
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_approx_d_sin.png", dpi=150)


    d2y  = -sin.(x)
    d2yh = d2ϕ * yi
    err = norm(d2y - d2yh) / norm(d2y) * 100

    println("diff2 err (percent) = ", err)

    plt.clf()
    plt.plot(x, d2y, label="exact")
    plt.plot(x, d2yh, label="approx")
    plt.grid()
    plt.legend()
    plt.savefig("TEMP_approx_d2_sin.png", dpi=150)

end

test_main()