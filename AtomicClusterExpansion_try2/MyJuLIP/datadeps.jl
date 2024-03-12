
module Deps 

using DataDeps

function fetch_test_pots() 
   register(DataDep(
         "MyJuLIP_testpots",
         "A few EAM potentials for testing",
         "https://www.dropbox.com/s/leub1c9ft1mm9fg/MyJuLIP_data.zip?dl=1",
         post_fetch_method = file -> run(`unzip $file`)
         ))
   return joinpath(datadep"MyJuLIP_testpots", "MyJuLIP_data") * "/"
end 


end