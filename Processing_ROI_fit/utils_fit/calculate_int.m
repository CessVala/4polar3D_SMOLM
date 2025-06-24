function I_fit = calculate_int(K, rho, eta, delta)
    M = M_iso(rho, eta, delta);
    I_fit = K*M;
    I_fit = I_fit./sum(I_fit);
end