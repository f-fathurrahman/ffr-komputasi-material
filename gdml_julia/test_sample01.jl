using Serialization: deserialize
using Statistics: quantile
using StatsBase: sample

E_all = deserialize("E.dat")

Nsamples = 200

p25 = quantile(E_all, 0.25)
p75 = quantile(E_all, 0.75)

h = 2 * (p75 - p25) / cbrt(Nsamples)
print("h = ", h)

Emax = maximum(E_all)
Emin = minimum(E_all)

if h > 0.0
    Nbins = Int64( ceil( (Emax - Emin)/h ) )
else
    Nbins = 1    
end

# Limit number of bins to half of requested subset size.
Nbins = min(Nbins, round(Int64, Nsamples/2))

bins = collect( range(Emin, Emax, Nbins+1)[1:Nbins] )
idxs = searchsortedlast.(Ref(bins), E_all)

# Get unique values
uniq_all = unique(idxs)

# Count occurrences of each unique value
cnts_all = [count(==(u), idxs) for u in uniq_all]

# limit reduced_cnts to what is available in cnts_all
reduced_cnts = Int64.( ceil.( cnts_all/sum(cnts_all) * Nsamples ) )

reduced_cnts = min.(reduced_cnts, cnts_all)

reduced_cnts_delta = Nsamples - sum(reduced_cnts)

# ??? only satisfied for reduced_cnts_delta == 0
#while abs(reduced_cnts_delta) > 0

    idxsg1 = reduced_cnts .> 1
    max_bin_reduction = min(reduced_cnts[idxsg1]) - 1

#end