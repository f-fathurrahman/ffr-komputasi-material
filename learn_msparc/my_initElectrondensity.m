function S = my_initElectrondensity(S)
if (S.ForceCount == 1)
    % use sum of the atomic charge densities as rho_guess
    S.rho = S.rho_at;
else
    % perform charge extrapolation
    if ((S.RelaxFlag == 1 || S.MDFlag) && S.ForceCount >= 4)
        S.rho(:,1) = S.rho_at(:,1) + S.delta_rho_in_tp1;
        S.rho(S.rho(:,1) < 0,1) = S.xc_rhotol;
        % update spin up/down densities
        if S.spin_typ > 0
            S.rho(:,2) = (S.rho(:,1) + S.mag(:,1)) * 0.5;
            S.rho(:,3) = (S.rho(:,1) - S.mag(:,1)) * 0.5; 
        end
    end
    % need to scale the density
    S.rho = S.rho * abs(S.NegCharge/dot(S.W,S.rho(:,1)));
end

end