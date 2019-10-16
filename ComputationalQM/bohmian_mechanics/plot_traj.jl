import PyPlot
const plt = PyPlot

using DelimitedFiles

function main()
    data = readdlm("TEMP_trajectory.dat")

    @views t = data[:,1]
    Np = size(data,2) - 1
    Nt = size(data,1)

    Nt2 = round(Int,Nt/2)

    println("Np = ", Np)
    plt.clf()
    for i in 1:Np
        plt.plot(t, data[:,i+1], color="black")
    end
    plt.savefig("TEMP_traj_x.png", dpi=150)
end

main()