diamond = Crystal(
    Lattice(2.527Å, 2.527Å, 2.527Å, 60, 60, 60),
    UnitCell(:C, [0.0, 0.0, 0.0], [1/4, 1/4, 1/4])
);

fig3d = Figure()
axds = Axis3(fig3d, xlabel="x/a₀", ylabel="y/a₀", zlabel="z/a₀", xgridvisible=false, ygridvisible=false, zgridvisible=false, aspect = :data)
axds_single = Axis3(fig3d, xlabel="x/a₀", ylabel="y/a₀", zlabel="z/a₀", xgridvisible=false, ygridvisible=false, zgridvisible=false, aspect = :data)

plotcrystal!(axds, diamond; ncells=2, showcell=false, showbonds=true)
plotcrystal!(axds_single, diamond; ncells=1, showcell=true, showbonds=true)
fig3d[1,1] = axds
fig3d[1,2] = axds_single

current_figure()


function kpath(kpoints, dk)
    vertices = reverse([point.second for point ∈ kpoints])
    labels = [point.first for point ∈ kpoints]
    path = [last(vertices)]
    plength = zeros(typeof(dk),1)
    idxs = [1]

    while length(vertices) >= 2
        pop!(path)
        v1 = pop!(vertices)
        v2 = last(vertices)
        dir = v2 .- v1
        dirm = norm(dir)
        segment = [v1 .+ dir .* dd 
            for dd ∈ range(start = 0, stop = 1, step = around(dk/dirm))]
        path = append!(path,segment)
        idxs = push!(idxs, last(idxs) + length(segment) - 1)
    end
    
    plength = append!(plength,
        Iterators.accumulate(
            +,[norm(v2 .- v1) for (v1,v2) ∈ zip(path[1:end-1], path[2:end])]
            ))
    points = [lab => plength[i] for (lab,i) ∈ zip(labels, idxs)]
    (path=path, plength=plength, ppoints=points)
end;



function elH(k, n, crystal::Crystal)
    G = crystal.lattice.G
    map(sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm)) do Gn
        q = G' * collect(k .- Gn)
        1/2*norm(q)^2
    end
end;


crystal1D = Crystal(
    Lattice(1.0),
    UnitCell(0.0)
)

kpath1D = -1.5:0.01:1.5

elH1D(k) = elH(k, 1, crystal1D)

es1D = elH1D.(kpath1D);

fig_1D = Figure()
ax_1D1 = Axis(fig_1D, title="Extended zone scheme")
ax_1D2 = Axis(fig_1D, title="First Brillouin zone")

for n ∈ 1:length(es1D[1])
    lines!(ax_1D1, kpath1D, [e[n] for e ∈ es1D])
    lines!(ax_1D2, kpath1D, [e[n] for e ∈ es1D])
end

hideydecorations!(ax_1D1)
hideydecorations!(ax_1D2)
ax_1D1.xticks = ([-1/2,1/2], ["-π/a","π/a"])
ax_1D2.xticks = ([-1/2,1/2], ["-π/a","π/a"])
ylims!(ax_1D1, (0,30))
ylims!(ax_1D2, (0,30))
xlims!(ax_1D2, (-1/2,1/2))

fig_1D[1,2] = ax_1D1
fig_1D[1,1] = ax_1D2

current_figure()

Al = Crystal(
    Lattice(2.856Å,2.856Å,2.856Å,60,60,60),
    UnitCell(:Al, [0.0,0.0,0.0])
);

Alks = kpath([
        :Γ => [0,0,0],
        :X => [1/2,0,1/2],
        :W => [1/2,1/4,3/4],
        :K => [3/8,3/8,3/4],
        :Γ => [0,0,0],
        :L => [1/2,1/2,1/2],
        :U => [5/8,1/4,5/8],
        :W => [1/2,1/4,3/4],
        :L => [1/2,1/2,1/2],
        :K => [3/8,3/8,3/4]
], 0.01);


alH_el(k) = elH(k, 2, Al)

es_al_el = alH_el.(Alks.path);


fig_Al_el = Figure()
ax_Al_el = Axis(fig_Al_el)
ax_Al_el.xticks = ([p.second for p ∈ Alks.ppoints], 
                 [string(p.first) for p ∈ Alks.ppoints])

ylims!(ax_Al_el, (0,1))
xlims!(ax_Al_el, (0, Alks.plength[end]))
hideydecorations!(ax_Al_el)

for n ∈ 1:length(es_al_el[1])
    lines!(Alks.plength, [e[n] for e ∈ es_al_el])
end

fig_Al_el[1,1] = ax_Al_el
current_figure()

#
# Perturbation theory
#

function nfH(k, n::Integer, V::Function, crystal::Crystal)
    G = crystal.lattice.G
    e = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm) |> collect
    k = G' * k
    Gs = (G' * collect(g) for g ∈ e)
    H = [V(j - i) for i ∈ Gs, j ∈ Gs]
    H .+= diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end

function nfH(k, e::Vector, V::Matrix, crystal::Crystal)
    G = crystal.lattice.G
    k = G' * k
    Gs = (G' * g for g ∈ e)
    V .+ diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end;

V1d = [0 1 1; 
       1 0 1;
       1 1 0]
nfH1D(k) = nfH(k, [-1,0,1], V1d, crystal1D)

nf1Des = nfH1D.(kpath1D) .|> eigvals;

fig_1D_nf = Figure()
ax_1D_nf = Axis(fig_1D_nf)

for n ∈ 1:length(es1D[1])
    lines!(ax_1D_nf, kpath1D, [e[n] for e ∈ es1D], color=:gray)
end

for n ∈ 1:length(nf1Des[1])
    lines!(ax_1D_nf, kpath1D, [e[n] for e ∈ nf1Des])
end

hideydecorations!(ax_1D_nf)
ax_1D_nf.xticks = ([-1/2,1/2], ["-π/a","π/a"])
ylims!(ax_1D_nf, (0,30))
xlims!(ax_1D_nf, (-1/2,1/2))

fig_1D_nf[1,1] = ax_1D_nf

current_figure()

# Al, screened potential
V(g, Q, q) = ifelse(norm(g) ≈ 0, 0, 4π * Q/(norm(g)^2 .+ q^2));
AlV(k) = V(k, 3, 10)
AlH(k) = nfH(k, 2, AlV, Al);
