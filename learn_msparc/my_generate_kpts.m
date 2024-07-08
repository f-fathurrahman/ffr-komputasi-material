function [S] = my_generate_kpts(S)
    nkpt = S.nkpt;
    if (S.BCx == 1 && nkpt(1) > 1)
        error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (x)');
    end
    if (S.BCy == 1 && nkpt(2) > 1)
        error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (y)');
    end
    if (S.BCz == 1 && nkpt(3) > 1)
        error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (z)');
    end

    % Monkhorst-pack grid for Brillouin zone sampling
%     MPG_typ1 = @(nkpt) (2*(1:nkpt) - nkpt - 1)/2; % MP grid points for infinite group order
    MPG_typ1 = @(nkpt) (-floor((nkpt - 1)/2):(-floor((nkpt - 1)/2)+nkpt-1));
    MPG_typ2 = @(nkpt) (0:nkpt-1); % MP grid points for finite group order

    if S.cell_typ < 3
        kptgrid_x = (1/nkpt(1)) * MPG_typ1(nkpt(1));
        kptgrid_y = (1/nkpt(2)) * MPG_typ1(nkpt(2));
        kptgrid_z = (1/nkpt(3)) * MPG_typ1(nkpt(3));
        sumx = 0;
        sumy = 0; 
        sumz = 0;
        % shift kpoint grid 
        kptgrid_x = kptgrid_x + S.kptshift(1) * (1/nkpt(1));
        kptgrid_y = kptgrid_y + S.kptshift(2) * (1/nkpt(2));
        kptgrid_z = kptgrid_z + S.kptshift(3) * (1/nkpt(3));
    
        % map k-points back to BZ
        temp_epsilon = eps; % include the right boundary k-points instead of left
        kptgrid_x = mod(kptgrid_x + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
        kptgrid_y = mod(kptgrid_y + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
        kptgrid_z = mod(kptgrid_z + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
    elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
        kptgrid_x = (1/nkpt(1)) * MPG_typ1(nkpt(1));
        kptgrid_y = (1/nkpt(2)) * MPG_typ2(nkpt(2));
        kptgrid_z = (1/nkpt(3)) * MPG_typ1(nkpt(3));
        sumx = 0;
        sumy = nkpt(2); 
        sumz = 0;
    end    
    
    % Scale kpoints
    kptgrid_x = (2*pi/S.L1) * kptgrid_x;
    kptgrid_y = (2*pi/S.L2) * kptgrid_y;
    kptgrid_z = (2*pi/S.L3) * kptgrid_z;

    [kptgrid_X, kptgrid_Y, kptgrid_Z] = ndgrid(kptgrid_x,kptgrid_y,kptgrid_z);
    kptgrid = [reshape(kptgrid_X,[],1),reshape(kptgrid_Y,[],1),reshape(kptgrid_Z,[],1)];
    disp(' reduced kpoint grid before symmetry:');
    disp(kptgrid*diag([S.L1/2/pi,S.L2/2/pi,S.L3/2/pi]));
    
    tnkpt = prod(nkpt);
    wkpt = ones(tnkpt,1)/tnkpt;% weights for k-points
    TOL = 1e-8;
    % Time-Reversal Symmetry to reduce k-points
    if S.TimeRevSym == 1
        Ikpt = zeros(tnkpt,1);
        Ikpt_rev = zeros(tnkpt,1);
        for ii = 1:tnkpt
            for jj = ii+1:tnkpt
                if (abs(kptgrid(ii,1) + kptgrid(jj,1) - sumx) < TOL) && (abs(kptgrid(ii,2) + kptgrid(jj,2) - sumy) < TOL) && (abs(kptgrid(ii,3) + kptgrid(jj,3) - sumz) < TOL)
                    Ikpt(ii) = 1;
                    Ikpt_rev(jj) = 1;
                end
            end
        end
        Ikpt = Ikpt>0.5;
        Ikpt_rev = Ikpt_rev>0.5;
        wkpt(Ikpt_rev) = 2*wkpt(Ikpt_rev);
        kptgrid = kptgrid(~Ikpt,:);
        wkpt = wkpt(~Ikpt);
        tnkpt = size(wkpt,1);
    end

    disp(' reduced kpoint grid after symmetry:');    
    disp(kptgrid*diag([S.L1/2/pi,S.L2/2/pi,S.L3/2/pi]));
    % Store into the structure
    S.kptgrid = kptgrid;
    S.tnkpt   = tnkpt;
    S.wkpt    = wkpt;
    
    % Generate kpoints grid for fock exchange
    if S.usefock == 1
        S.isgamma = 0;
        if tnkpt == 1 && sum(kptgrid == [0,0,0])==3
            S.isgamma = 1;
        end
        
        % Use part of full k-point grid
        if S.exx_downsampling(1) == 0
            kptgrid_x_hf = 0;
            if sum(find(ismembertol(kptgrid_x,0,1e-8))) == 0
                error("Gamma point is not one of the k-vectors. Please use positive EXX_DOWNSAMPLING or change k-point grid in first direction.");
            end
        else
            range = S.exx_downsampling(1):S.exx_downsampling(1):nkpt(1);
            kptgrid_x_hf = kptgrid_x(range);
        end
        
        if S.exx_downsampling(2) == 0
            kptgrid_y_hf = 0;
            if sum(find(ismembertol(kptgrid_y,0,1e-8))) == 0
                error("Gamma point is not one of the k-vectors. Please use positive EXX_DOWNSAMPLING or change k-point grid in second direction.");
            end
        else
            range = S.exx_downsampling(2):S.exx_downsampling(2):nkpt(2);
            kptgrid_y_hf = kptgrid_y(range);
        end
        
        if S.exx_downsampling(3) == 0
            kptgrid_z_hf = 0;
            if sum(find(ismembertol(kptgrid_z,0,1e-8))) == 0
                error("Gamma point is not one of the k-vectors. Please use positive EXX_DOWNSAMPLING or change k-point grid in third direction.");
            end
        else
            range = S.exx_downsampling(3):S.exx_downsampling(3):nkpt(3);
            kptgrid_z_hf = kptgrid_z(range);
        end
        
        [kptgrid_X_HF, kptgrid_Y_HF, kptgrid_Z_HF] = ndgrid(kptgrid_x_hf,kptgrid_y_hf,kptgrid_z_hf);
        kptgrid_HF = [reshape(kptgrid_X_HF,[],1),reshape(kptgrid_Y_HF,[],1),reshape(kptgrid_Z_HF,[],1)];
        disp(' reduced kpoint grid for Fock Exchange operator.');
        disp(kptgrid_HF*diag([S.L1/2/pi,S.L2/2/pi,S.L3/2/pi]));
        
        S.kptgridhf = kptgrid_HF;
        S.tnkpthf   = length(kptgrid_x_hf)*length(kptgrid_y_hf)*length(kptgrid_z_hf);
        S.wkpthf    = ones(S.tnkpthf,1)/S.tnkpthf;
        S.nkpthf = [length(kptgrid_x_hf), length(kptgrid_y_hf), length(kptgrid_z_hf)];
    end
end