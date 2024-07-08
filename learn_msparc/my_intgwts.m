function [W] = my_intgwts(Nx,Ny,Nz,BCx,BCy,BCz,xin,S)
    if S.cell_typ == 1 || S.cell_typ == 2
        W_x = ones(Nx,1)*S.dx;
        W_x(1) = W_x(1) * (1-BCx*0.5);
        W_x(Nx) = W_x(Nx) * (1-BCx*0.5);

        W_y = ones(Ny,1)*S.dy;
        W_y(1) = W_y(1) * (1-BCy*0.5);
        W_y(Ny) = W_y(Ny) * (1-BCy*0.5);

        W_z = ones(Nz,1)*S.dz;
        W_z(1) = W_z(1) * (1-BCz*0.5);
        W_z(Nz) = W_z(Nz) * (1-BCz*0.5);

        W = kron(W_z,kron(W_y,W_x)) * S.Jacb;
    elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
        W = IntgWts_cychel(Nx,Ny,Nz,BCx,BCy,BCz,xin,S);
        % Not called ??
    end
end