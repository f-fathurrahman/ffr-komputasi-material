# Example arrays
x = [0.2, 6.4, 3.0, 1.6]
bins = [0.0, 1.0, 2.5, 4.0, 10.0]

# Using searchsortedlast with broadcasting to digitize
indices = searchsortedlast.(Ref(bins), x)


# Unique return counts

# Example array
idxs = [1, 2, 2, 3, 4, 4, 4, 5]

# Get unique values
unique_vals = unique(idxs)

# Count occurrences of each unique value
counts = [count(==(u), idxs) for u in unique_vals]

# Combine unique values and their counts into a dictionary
unique_counts = Dict(unique_vals .=> counts)
