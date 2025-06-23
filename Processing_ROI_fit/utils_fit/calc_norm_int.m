function I_norm = calc_norm_int(K, M)
    I_norm = K*M;
    I_norm = I_norm./sum(I_norm,1);
end