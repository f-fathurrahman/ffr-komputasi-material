include("inc_crystal.jl")

hBN = Crystal(
    Lattice(2.512Å, 2.512Å, 120),
    UnitCell([:B, :N], [2/3, 1/3], [1/3, 2/3])
)

fig2ds = Figure()
axhbns = Axis(fig2ds;
    title = "hBN",
    xlabel = "x/a₀",
    ylabel = "y/a₀",
    xgridvisible = false,
    ygridvisible = false,
    aspect = DataAspect()
)

plotcrystal!(axhbns, hBN, ncells=3)
fig2ds[1,1] = axhbns

current_figure()
