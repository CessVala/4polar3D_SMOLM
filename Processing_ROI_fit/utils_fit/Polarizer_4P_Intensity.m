function I_pol = Polarizer_4P_Intensity(n_sample, n_glass, NA, bfp_factor, alpha, z_dipole, lambda, SAF_flag, NAred_ch,varargin)
    
    if nargin == 10
        A = varargin{1};
    else
        A = ones(4,1);
    end

    n = n_glass/n_sample;                      % numerical aperture
    
    if SAF_flag
        theta_max = asin(NA./n_glass);
        theta_red = asin(bfp_factor.*NA./n_glass);
    else
        theta_max = asin(n_sample./n_glass);
        theta_red = asin(bfp_factor.*n_sample./n_glass);
    end

    % General parameters
    
    cos2        = @(p1,t1) sqrt(1-(n_glass/n_sample.*sin(t1)).^2);
    sin_2phi1   = @(p1,t1) sin(2*p1);
    rho         = @(p1,t1) sin(t1);
    
    t_p = @(p1,t1) 2*n_sample*cos2(p1,t1)./(n_sample.*cos(t1)+n_glass.*cos2(p1,t1));     % Fresnel transmission coefficient for p-polar
    t_s = @(p1,t1) 2*n_sample*cos2(p1,t1)./(n_sample.*cos2(p1,t1)+n_glass.*cos(t1));     % Fresnel transmission coefficient for s-polar
    
    psy_depth   = @(p1,t1) 2*pi*z_dipole*n_sample/lambda*cos2(p1,t1);

    % Electric fields to integrate
    
    E_x_mux     = @(p1,t1)  n   .*  ( t_s(p1,t1) .*cos(t1)./cos2(p1,t1).*sin(p1).^2  + t_p(p1,t1) .*sqrt( 1-rho(p1,t1).^2 ).*cos(p1).^2 ) .* exp( 1i*psy_depth(p1,t1) );
    E_x_muy     = @(p1,t1) -n/2 .*  ( t_s(p1,t1) .*cos(t1)./cos2(p1,t1)              - t_p(p1,t1) .*sqrt( 1-rho(p1,t1).^2 ) )             .* exp( 1i*psy_depth(p1,t1) ).*sin_2phi1(p1,t1);
    E_x_muz     = @(p1,t1) -n^2 .*  ( t_p(p1,t1) .*cos(t1)./cos2(p1,t1).*rho(p1,t1).*cos(p1) )                                            .* exp( 1i*psy_depth(p1,t1) );
    
    E_y_mux     = @(p1,t1) -n/2 .*  ( t_s(p1,t1) .*cos(t1)./cos2(p1,t1)              - t_p(p1,t1) .*sqrt( 1-rho(p1,t1).^2 ) )             .* exp( 1i*psy_depth(p1,t1) ).*sin_2phi1(p1,t1);
    E_y_muy     = @(p1,t1)  n   .*  ( t_s(p1,t1) .*cos(t1)./cos2(p1,t1).*cos(p1).^2  + t_p(p1,t1) .*sin(p1).^2.*sqrt(1-rho(p1,t1).^2) )   .* exp(1i*psy_depth(p1,t1));
    E_y_muz     = @(p1,t1) -n^2 .*  ( t_p(p1,t1) .*cos(t1)./cos2(p1,t1).*rho(p1,t1).*sin(p1) )                                            .* exp( 1i*psy_depth(p1,t1) );
       
    E_45_mux    = @(p1,t1) ( E_x_mux(p1,t1) + E_y_mux(p1,t1) )/sqrt(2);
    E_45_muy    = @(p1,t1) ( E_x_muy(p1,t1) + E_y_muy(p1,t1) )/sqrt(2);
    E_45_muz    = @(p1,t1) ( E_x_muz(p1,t1) + E_y_muz(p1,t1) )/sqrt(2);
    
    E_135_mux   = @(p1,t1) (-E_x_mux(p1,t1) + E_y_mux(p1,t1) )/sqrt(2);
    E_135_muy   = @(p1,t1) (-E_x_muy(p1,t1) + E_y_muy(p1,t1) )/sqrt(2);
    E_135_muz   = @(p1,t1) (-E_x_muz(p1,t1) + E_y_muz(p1,t1) )/sqrt(2);

    % Factors of the matrix
    m000   = @(p1,t1) abs(E_x_mux(p1,t1)).^2 .* sin(t1) + abs(E_x_muy(p1,t1)).^2 .* sin(t1) + abs(E_x_muz(p1,t1)).^2  .* sin(t1);
    m090   = @(p1,t1) abs(E_y_mux(p1,t1)).^2 .* sin(t1) + abs(E_y_muy(p1,t1)).^2 .* sin(t1) + abs(E_y_muz(p1,t1)).^2  .* sin(t1);
    m045   = @(p1,t1) abs(E_45_mux(p1,t1)).^2 .* sin(t1) + abs(E_45_muy(p1,t1)).^2 .* sin(t1) + abs(E_45_muz(p1,t1)).^2  .* sin(t1);
    m135   = @(p1,t1) abs(E_135_mux(p1,t1)).^2 .* sin(t1) + abs(E_135_muy(p1,t1)).^2 .* sin(t1) + abs(E_135_muz(p1,t1)).^2  .* sin(t1);

   
      
    %% theoretical K matrix from individual dipoles directions

    if NAred_ch == 1
        % intensity integrals high NA channels
       
        M045 = integral2(m045,0,2*pi,0,theta_max,'method','iterated');        
        M135 = integral2(m135,0,2*pi,0,theta_max,'method','iterated');
    
        % intensity integrals low NA channels
        M000 = integral2(m000,0,2*pi,0,theta_red,'method','iterated');        
        M090 = integral2(m090,0,2*pi,0,theta_red,'method','iterated');

    elseif NAred_ch == 2
        % intensity integrals high NA channels
        M000 = integral2(m000,0,2*pi,0,theta_max,'method','iterated');        
        M090 = integral2(m090,0,2*pi,0,theta_max,'method','iterated');

        % intensity integrals low NA channels
        M045 = integral2(m045,0,2*pi,0,theta_red,'method','iterated');        
        M135 = integral2(m135,0,2*pi,0,theta_red,'method','iterated');

    else
        errordlg('Wrong choice for NAred_ch! Choose either 1 (0/90 reduced NA) or 2 (45/135 reduced NA)');
        return
    end
    % REMARK: This is not the same as a K matrix
    Pol_Mat     = [M000;
                   M090;
                   M045; 
                   M135];
    
    I_pol       = Pol_Mat .*(1 - A.*sin(alpha+[0; pi/2; -pi/4; +pi/4]).^2);
    I_pol       = I_pol./sum(I_pol,1);
end