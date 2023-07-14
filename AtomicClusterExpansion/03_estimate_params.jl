solver = ACEfit.BLR() 
my_acefit!(model, train; solver=solver, data_keys...);

# To see the training errors we can use 

@info("Training Errors")
my_linear_errors(train, model; data_keys...);