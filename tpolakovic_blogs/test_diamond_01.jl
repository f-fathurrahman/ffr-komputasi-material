include("inc_crystal.jl")

diamond = Crystal(
    Lattice(2.527Å, 2.527Å, 2.527Å, 60, 60, 60),
    UnitCell(:C, [0.0, 0.0, 0.0], [1/4, 1/4, 1/4])
);

fig3d = Figure()
axds = Axis3(fig3d,
    xlabel="x/a₀",
    ylabel="y/a₀",
    zlabel="z/a₀",
    xgridvisible=false,
    ygridvisible=false,
    zgridvisible=false,
    aspect=:data)

axds_single = Axis3(fig3d, xlabel="x/a₀", ylabel="y/a₀", zlabel="z/a₀", xgridvisible=false, ygridvisible=false, zgridvisible=false, aspect = :data)

plotcrystal!(axds, diamond; ncells=2, showcell=false, showbonds=true)
plotcrystal!(axds_single, diamond; ncells=1, showcell=true, showbonds=true)

fig3d[1,1] = axds
fig3d[1,2] = axds_single

current_figure()
