function M_out = M_iso(rho, eta, delta)
   
    l_iso   =   lambda_iso(delta);
    
    M_out   =   [(1 - 3.*l_iso) .* cos(rho).^2 .* sin(eta).^2 + l_iso; ...
                 (1 - 3.*l_iso) .* sin(rho).^2 .* sin(eta).^2 + l_iso; ...
                 (1 - 3.*l_iso) .* cos(eta).^2 + l_iso; ...
                 (1 - 3.*l_iso) .* sin(eta).^2 .* cos(rho) .* sin(rho)];
end