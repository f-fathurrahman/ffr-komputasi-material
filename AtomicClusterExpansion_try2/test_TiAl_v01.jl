using MyACEpotentials

data = read_extxyz("datasets/TiAl_tutorial.xyz")
# we don't need other data, only train

model = acemodel(
    elements = [:Ti, :Al],
    order = 3, totaldegree = 12,
    Eref = [:Ti => -1586.0195, :Al => -105.5954]
)

acefit!(model, data);
errors = MyACEpotentials.linear_errors(data, model)
