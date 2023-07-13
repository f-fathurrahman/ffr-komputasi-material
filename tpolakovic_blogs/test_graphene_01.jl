include("inc_crystal.jl")

graphene = Crystal(
    Lattice(2.468Å, 2.468Å, 120),
    UnitCell(:C, [2/3, 1/3], [1/3, 2/3])
)

# graphene.lattice.R[:,1] ⋅ graphene.lattice.R[:,2] != 0

fig2ds = Figure()
axgs = Axis(fig2ds;
    title = "graphene",
    xlabel = "x/a₀",
    ylabel = "y/a₀",
    xgridvisible = false,
    ygridvisible = false,
    aspect = DataAspect()
)

plotcrystal!(axgs, graphene, ncells=2)
fig2ds[1,1] = axgs
current_figure()
