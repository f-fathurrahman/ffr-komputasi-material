# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Setup

# %%
include("setup_works.jl")

# %%
using MyGaussianBasis

# %% [markdown]
# HCl molecule

# %% [markdown]
# Construct basis set for a given molecule (provided as string):

# %%
bs_HCl = BasisSet("sto-3g", """
H         0.00      0.00     0.00
Cl        0.76      0.00     0.00""")

# %%
typeof(bs_HCl)

# %%
fieldnames(typeof(bs_HCl))

# %%
length(bs_HCl.basis) # one for each atoms

# %%
bs_HCl.basis[1] |> typeof

# %% [markdown]
# There is only one basis functio for H atom:

# %%
length(bs_HCl.basis[1])

# %% [markdown]
# There are five basis functions for Cl:

# %%
length(bs_HCl.basis[2])

# %%
for (i,b) in enumerate(bs_HCl.basis[1])
    println("\nBasis number: $i")
    display(b)
end

# %%
for (i,b) in enumerate(bs_HCl.basis[2])
    println("\nBasis number: $i")
    display(b)
end

# %% [markdown]
# ## Overlap matrix elements

# %%
𝒪 = overlap(bs_HCl)

# %%
MyGaussianBasis |> names
