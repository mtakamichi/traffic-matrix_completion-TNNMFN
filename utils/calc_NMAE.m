% Calc. Normalized Mean Absolute Error

function NMAE =calc_NMAE(T_org,T_mask,T_est)
    NMAE = sum(abs(T_org(:).*~T_mask(:)-T_est(:).*~T_mask(:)))/ sum(abs(T_org(:).*~T_mask(:)));
end