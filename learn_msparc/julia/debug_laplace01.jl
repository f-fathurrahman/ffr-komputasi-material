# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
Lx = 13.0;
Nx = 27; # Number of discretization points

# %%
FDn = 6 # finite difference order
# using FDn*2 + 1 points

# %%
dx = Lx/(Nx-1);

# %%
# Finite difference weights of the second derivative
w2 = zeros(Float64, FDn+1)
for k in 1:FDn
    w2[k+1] = (2*(-1)^(k+1))*(factorial(FDn)^2)/(k*k*factorial(FDn-k)*factorial(FDn+k))
    w2[1] = w2[1] - 2*(1/(k*k))
end

# %%
w2

# %%
# Finite difference weights of the first derivative
w1 = zeros(Float64, FDn+1) # should be FDn+1 not FDn
for k in 1:FDn
    w1[k+1] = ((-1)^(k+1))*(factorial(FDn)^2)/(k*factorial(FDn-k)*factorial(FDn+k))
end

# %% [markdown]
# 1st index of w1 is not used?

# %%
w1

# %%
n0 = FDn;

# %%
# Initial number of non-zeros: including ghost nodes
tmp_nnzCount = (2 * n0 + 1) * Nx;

# %%
tmp_nnzCount

# %%
# Row and column indices and the corresponding non-zero values
# used to generate sparse matrix DL11 s.t. DL11(I(k),II(k)) = V(k)
tmp_idx_row = zeros(Int64, tmp_nnzCount)
tmp_idx_col = zeros(Int64, tmp_nnzCount)
tmp_vv = zeros(Float64, tmp_nnzCount)
row_cnt = 1
cnt = 1
coef_dxx = 1/dx^2
# Find non-zero entries that use forward difference
for ii in 1:Nx
    # diagonal element
    tmp_idx_row[cnt] = row_cnt
    tmp_idx_col[cnt] = ii
    tmp_vv[cnt] = w2[1]*coef_dxx
    cnt += 1 # increment
    # off-diagonal elements
    for q in 1:n0
        # ii + q
        tmp_idx_row[cnt] = row_cnt
        tmp_idx_col[cnt] = ii + q
        tmp_vv[cnt] = w2[q+1]*coef_dxx
        cnt += 1 # increment
        # ii - q
        tmp_idx_row[cnt] = row_cnt
        tmp_idx_col[cnt] = ii - q
        tmp_vv[cnt] = w2[q+1]*coef_dxx
        cnt += 1 # increment
    end
    row_cnt += 1
end

# %%
BCx = :DIRICHLET_BC
if BCx == :DIRICHLET_BC
    # Removing outside domain entries
    idx_inside = (tmp_idx_col .>= 1) .&& (tmp_idx_col .<= Nx)
    idx_row = tmp_idx_row[idx_inside]
    idx_col = tmp_idx_col[idx_inside]
    vv = tmp_vv[idx_inside]
end;

# %%
AA = zeros(Float64, Nx, Nx)
nnzCount = length(I_11)
for ip in 1:nnzCount
    i = idx_row[ip]
    j = idx_col[ip]
    Aij = vv[ip]
    AA[i,j] = Aij
end

# %%
using Plots: heatmap, theme

# %%
using PlotThemes

# %%
theme(:dark)

# %%
heatmap(AA, yflip=true, aspect_ratio=:equal)

# %%
AA[1,2], AA[2,1]
