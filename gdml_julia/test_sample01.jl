using Serialization: deserialize
using Statistics: quantile
using StatsBase: sample, ProbabilityWeights

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
while abs(reduced_cnts_delta) > 0

    idxsg1 = reduced_cnts .> 1
    max_bin_reduction = minimum(reduced_cnts[idxsg1]) - 1

    probs = (reduced_cnts .- 1) / sum(reduced_cnts .- 1)
    NbinsAdditional = min(max_bin_reduction, abs.(reduced_cnts_delta))

    outstanding = sample(uniq_all, ProbabilityWeights(probs),
        NbinsAdditional, replace=true)

    uniq_outstanding = unique(outstanding)
    cnts_outstanding = [count(==(u), outstanding) for u in uniq_outstanding]

    # Bucket IDs to Idxs.
    outstanding_bucket_idx = findall(in.(uniq_outstanding, uniq_all))
    println("outstanding_bucket_idx = ", outstanding_bucket_idx)

    reduced_cnts[outstanding_bucket_idx] += (sign(reduced_cnts_delta) .* cnts_outstanding)
    print("Updated reduced_cnts")

    reduced_cnts_delta = Nsamples - sum(reduced_cnts)

end

idxs_train = Vector{Int64}()
for (uniq_idx, bin_cnt) in zip(uniq_all, reduced_cnts)
    println("uniq_idx = ", uniq_idx, " bin_cnt = ", bin_cnt)
    idx_in_bin_all = findall(idxs .== uniq_idx)
    append!(idxs_train, sample(idx_in_bin_all, bin_cnt, replace=false))
end

