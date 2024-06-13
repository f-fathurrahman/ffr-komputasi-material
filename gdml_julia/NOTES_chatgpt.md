In Julia, the equivalent of Numpy's `np.random.choice` can be achieved using the `sample` function from the `StatsBase.jl` package. This function allows you to perform random sampling with or without replacement and can also handle non-uniform probability distributions. Here's how you can use it:

```julia
using StatsBase

# Example array and probabilities
a = 1:5
p = [0.1, 0.2, 0.3, 0.2, 0.2]

# Sample with replacement
random_sample = sample(a, Weights(p), 3)

# Sample without replacement
random_sample_no_replace = sample(a, Weights(p), 3, replace=false)
```

In the above example, `random_sample` will give you a sample of 3 elements from the array `a` with the given probabilities `p`, with replacement. The `random_sample_no_replace` will do the same without replacement¹.

If you need to generate a non-uniform random sample based on a user-supplied probability list, you can use `ProbabilityWeights` instead of `Weights`:

```julia
# Example probability list
probability_list = [0.2, 0.2, 0.2, 0.2, 0.2]

# Generate random indices based on non-uniform random sampling
rand_idx = sample(1:300, ProbabilityWeights(probability_list))
```

This will generate random indices from 1 to 300 based on the supplied `probability_list`¹.

For cases where you want to use the `Distributions.jl` package, you can create a `Categorical` distribution and then use `rand` to sample from it:

```julia
using Distributions

# Create a categorical distribution
a = Categorical(p)

# Sample from the categorical distribution
rand_sample = rand(a)
```

This will return the index of the array with the associated probability¹.

Remember to add the `StatsBase` or `Distributions` package to your Julia environment if you haven't already done so. You can add them using Julia's package manager with the commands `import Pkg; Pkg.add("StatsBase")` and `import Pkg; Pkg.add("Distributions")`.

Source: Conversation with Copilot, 6/12/2024
(1) What is the equivalent of numpy.random.choice for generating a non .... https://discourse.julialang.org/t/what-is-the-equivalent-of-numpy-random-choice-for-generating-a-non-uniform-random-sample-based-on-a-user-supplied-probability/34461.
(2) numpy.random.choice — NumPy v1.26 Manual. https://numpy.org/doc/stable/reference/random/generated/numpy.random.choice.html.
(3) function to choose element at random from an array #3075 - GitHub. https://github.com/JuliaLang/julia/issues/3075.
(4) Julia: Construct an array from an index array and an array to choose .... https://stackoverflow.com/questions/77475843/julia-construct-an-array-from-an-index-array-and-an-array-to-choose-from-or-w.

