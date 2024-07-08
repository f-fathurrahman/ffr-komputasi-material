function S = my_calculate_rb(S)
% Starting and ending indices of b-region
if (S.cell_typ == 1 || S.cell_typ == 2)
    pos_atm_x = 0; % atom location in x-direction
    pos_atm_y = 0; % atom location in y-direction
    pos_atm_z = 0; % atom location in z-direction
    rb_up_x = (S.dx < 1.5) * (10+10*S.dx) + (S.dx >=1.5) * (20*S.dx-9.5);
    rb_up_y = (S.dy < 1.5) * (10+10*S.dy) + (S.dy >=1.5) * (20*S.dy-9.5);
    rb_up_z = (S.dz < 1.5) * (10+10*S.dz) + (S.dz >=1.5) * (20*S.dz-9.5);
    f_rby = @(y) y;
    
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
    pos_atm_x = S.xmax_at; % maximum R coordinate of any atom
    pos_atm_y = 0; % atom location in theta-direction
    pos_atm_z = 0; % atom location in z-direction
    rb_up_x = S.xvac; % Radial direction vacuum
    f_rby = @(y) acos(1 - y^2/(2*pos_atm_x^2));
    rb_up_y = f_rby(12); % Theta direction
    rb_up_z = 12; % z-direction
    
end
ii_s_temp = -ceil(rb_up_x/S.dx);
ii_e_temp = ceil(rb_up_x/S.dx);
jj_s_temp = -ceil(rb_up_y/S.dy);
jj_e_temp = ceil(rb_up_y/S.dy);
kk_s_temp = 0;
kk_e_temp = ceil(rb_up_z/S.dz);
xx_temp = pos_atm_x + (ii_s_temp-S.FDn:ii_e_temp+S.FDn)*S.dx;
yy_temp = pos_atm_y + (jj_s_temp-S.FDn:jj_e_temp+S.FDn)*S.dy;
zz_temp = pos_atm_z + (kk_s_temp-S.FDn:kk_e_temp+S.FDn)*S.dz;
[XX_3D_temp,YY_3D_temp,ZZ_3D_temp] = ndgrid(xx_temp,yy_temp,zz_temp);
Nx = (ii_e_temp-ii_s_temp)+1;
Ny = (jj_e_temp-jj_s_temp)+1;
Nz = (kk_e_temp-kk_s_temp)+1;
% Find distances
dd_temp = calculateDistance(XX_3D_temp,YY_3D_temp,ZZ_3D_temp,pos_atm_x,pos_atm_y,pos_atm_z,S);

% Find integration weights
W_temp = my_intgwts(Nx,Ny,Nz,1,1,1,xx_temp(S.FDn+1),S); % 1 - dirichlet BC on the boundary nodes
W_temp = reshape(W_temp,Nx,Ny,Nz);

% Find VJ and bJ
for ityp = 1:S.n_typ
    V_PS_temp = zeros(size(dd_temp));
    IsLargeThanRmax = dd_temp > S.Atm(ityp).r_grid_vloc(end);
    V_PS_temp(IsLargeThanRmax) = -S.Atm(ityp).Z;
    V_PS_temp(~IsLargeThanRmax) = interp1(S.Atm(ityp).r_grid_vloc, ...
        S.Atm(ityp).r_grid_vloc.*S.Atm(ityp).Vloc, dd_temp(~IsLargeThanRmax), 'spline');
    
    V_PS_temp = V_PS_temp./dd_temp;
    V_PS_temp(dd_temp<S.Atm(ityp).r_grid_vloc(2)) = S.Atm(ityp).Vloc(1);
    II_temp = 1+S.FDn : size(V_PS_temp,1)-S.FDn;
    JJ_temp = 1+S.FDn : size(V_PS_temp,2)-S.FDn;
    KK_temp = 1+S.FDn : size(V_PS_temp,3)-S.FDn;
    
    b_temp = pseudochargeDensity_atom(V_PS_temp,II_temp,JJ_temp,KK_temp,xx_temp(1),S);
    b_temp = -b_temp / (4*pi);
    err_rb = 100;
    count = 1;
    rb_x = S.Atm(ityp).rc;
    rb_y = f_rby(S.Atm(ityp).rc);
    rb_z = S.Atm(ityp).rc;
    rb_x = ceil(rb_x/S.dx-1e-12)*S.dx;
    rb_y = ceil(rb_y/S.dy-1e-12)*S.dy;
    rb_z = ceil(rb_z/S.dz-1e-12)*S.dz;
    fprintf(' Finding rb for %s ...\n',S.Atm(ityp).typ);
    while (err_rb > S.pseudocharge_tol && count <= 100 && rb_x <= rb_up_x && rb_y <= rb_up_y && rb_z <= rb_up_z )
        rb_x = rb_x + S.dx;
        rb_z = rb_z + S.dz;
        rb_y = f_rby(max(rb_x,rb_z));
        ii_rb = -1*ii_s_temp+S.FDn-floor(rb_x/S.dx)+1:-1*ii_s_temp+S.FDn+floor(rb_x/S.dx)+1;
        jj_rb = -1*jj_s_temp+S.FDn-floor(rb_y/S.dy)+1:-1*jj_s_temp+S.FDn+floor(rb_y/S.dy)+1;
        kk_rb = S.FDn+1:S.FDn+floor(rb_z/S.dz)+1; 
        err_rb = abs(sum(sum(sum(W_temp(ii_rb-S.FDn,jj_rb-S.FDn,kk_rb-S.FDn).*b_temp(ii_rb,jj_rb,kk_rb))))*2 + S.Atm(ityp).Z);
        fprintf(' rb = {%.3f %.3f %.3f}, int_b = %.15f, err_rb = %.3e\n',rb_x,rb_y,rb_z,2*sum(sum(sum(W_temp(ii_rb-S.FDn,jj_rb-S.FDn,kk_rb-S.FDn).*b_temp(ii_rb,jj_rb,kk_rb)))),err_rb);
        count = count + 1;
    end
    
    assert(rb_x<=rb_up_x && rb_y<=rb_up_y && rb_z<=rb_up_z,'Need to increase upper bound for rb!');
    S.Atm(ityp).rb_x = rb_x;
    S.Atm(ityp).rb_y = rb_y;
    S.Atm(ityp).rb_z = rb_z;
    % S.Atm(ityp).rb_x = ceil(rb_x/S.dx-1e-12)*S.dx; % + S.dx;
    % S.Atm(ityp).rb_y = ceil(rb_y/S.dy-1e-12)*S.dy; % + S.dy;
    % S.Atm(ityp).rb_z = ceil(rb_z/S.dz-1e-12)*S.dz; % + S.dz;
    fprintf(' rb = {%.3f %.3f %.3f}\n',S.Atm(ityp).rb_x,S.Atm(ityp).rb_y,S.Atm(ityp).rb_z);
end

end