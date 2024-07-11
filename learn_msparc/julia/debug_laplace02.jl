# Case of periodic, non-orthonal cell

include("create_D1_matrix.jl")
include("create_D2_matrix.jl")

BCx = :PERIODIC_BC

LatVecs = zeros(Float64, 3, 3) # !!! stored by rows !!!!
LatVecs[1,:] = [0.667124384994991, 0.741249316661101, 0.074124931666110]*6
LatVecs[2,:] = [0.074124931666110, 0.741249316661101, 0.667124384994991]*6
LatVecs[3,:] = [0.667124384994991, 0.074124931666110, 0.741249316661101]*6
# LatVecs scale: 6


# Check the cell typ (orthogonal or non-orthogonal)
if(abs(dot(S.lat_vec(1,:),S.lat_vec(2,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(2,:),S.lat_vec(3,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(3,:),S.lat_vec(1,:))) > S.temp_tol)
    S.cell_typ = 2;
end

S.lat_uvec(1,:) = S.lat_vec(1,:)/norm(S.lat_vec(1,:));
S.lat_uvec(2,:) = S.lat_vec(2,:)/norm(S.lat_vec(2,:));
S.lat_uvec(3,:) = S.lat_vec(3,:)/norm(S.lat_vec(3,:));

% Set up transformation matrices for non orthogonal cells
if S.cell_typ == 2
    % Jacobian
    S.Jacb = det(S.lat_uvec');
    assert(S.Jacb > 0.0,'Volume is negative!');

    % metric_T, Gradient and laplacian transformation matrices
    S.metric_T = S.lat_uvec * S.lat_uvec' ;
    S.metric_T(1,2) = 2*S.metric_T(1,2); 
    S.metric_T(2,3) = 2*S.metric_T(2,3); 
    S.metric_T(1,3) = 2*S.metric_T(1,3);
    S.grad_T = inv(S.lat_uvec') ;
    S.lapc_T = S.grad_T * S.grad_T' ;
    S.lapc_T(1,2) = 2*S.lapc_T(1,2); 
    S.lapc_T(2,3) = 2*S.lapc_T(2,3);
    S.lapc_T(1,3) = 2*S.lapc_T(1,3);

    count_prev = 0;
    count = 0;
    for ityp = 1:S.n_typ
        if(S.IsFrac(ityp) == 0)
            S.Atm(ityp).coords = transpose(S.grad_T * transpose(S.Atm(ityp).coords));
            count = count + S.Atm(ityp).n_atm_typ;
            S.Atoms(count_prev+1:count,:) = S.Atm(ityp).coords;
            count_prev = count;
        else
            count = count + S.Atm(ityp).n_atm_typ;
            count_prev = count;
        end
    end
end



Lx = 6.0
Nx = 15 # Number of discretization points
FDn = 6 # finite difference order
# using FDn*2 + 1 points

D2mat = create_D2_matrix_FD(FDn, Nx, :PERIODIC_BC)
D1mat = create_D1_matrix_FD(FDn, Nx, :PERIODIC_BC)

#=
using Plots: heatmap, theme
using PlotThemes
theme(:dark)
heatmap(D2mat, yflip=true, aspect_ratio=:equal)
=#

