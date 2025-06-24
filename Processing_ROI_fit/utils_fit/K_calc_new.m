function K = K_calc_new(bfp_factor, n_glass, n_sample, NA, SAF_flag, varargin)
    
    if nargin == 5
        lambda      = 0.600;        % wavelength (in um)
        z_dipole    = 0.0;
        red_ch_select = 1;
    elseif nargin == 6
        lambda      = 0.600;        % wavelength (in um)
        z_dipole    = varargin{1};
        red_ch_select = 1;
    elseif nargin == 7
        lambda      = varargin{2};  % wavelength (in um)
        z_dipole    = varargin{1};
        red_ch_select = 1;
    elseif nargin == 8
        lambda      = varargin{2};  % wavelength (in um)
        z_dipole    = varargin{1};
        red_ch_select = varargin{3};
    else
        errordlg('Wrong Inputs!!!')
        return
    end

    n = n_glass/n_sample;                      % numerical aperture
    
    if SAF_flag
        theta_max = asin(NA./n_glass);
        theta_red = asin(bfp_factor.*NA./n_glass);
    else
        theta_max = asin(n_sample./n_glass);
        theta_red = asin(bfp_factor.*n_sample./n_glass);
    end
    %%

    %%
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
    
    % E_x_mux     = @(p1,t1)  n   .*  ( ( cos(t1)./cos2(p1,t1) ) .* t_s(p1,t1) .* sin(p1).^2 + t_p(p1,t1) .* cos(p1).^2 .* sqrt(1-rho(p1,t1).^2)) .* exp(1i*psy_depth(p1,t1));
    % E_x_muy     = @(p1,t1) -n/2 .*  sin_2phi1(p1,t1).* ( (cos(t1)./cos2(p1,t1) ) .* t_s(p1,t1) - t_p(p1,t1) .* sqrt(1-rho(p1,t1).^2)) .* exp(1i*psy_depth(p1,t1));
    % E_x_muz     = @(p1,t1) -n^2 .*  (cos(t1)./cos2(p1,t1)).*t_p(p1,t1).*rho(p1,t1).*cos(p1).* exp(1i*psy_depth(p1,t1));
    % 
    % E_y_mux     = @(p1,t1) -n/2 .*  sin_2phi1(p1,t1) .* ( (cos(t1)./cos2(p1,t1)).*t_s(p1,t1) - t_p(p1,t1).*sqrt(1-rho(p1,t1).^2) ).* exp(1i*psy_depth(p1,t1));
    % E_y_muy     = @(p1,t1)  n   .*  ((cos(t1)./cos2(p1,t1)).*t_s(p1,t1).*cos(p1).^2 + t_p(p1,t1).*sin(p1).^2.*sqrt(1-rho(p1,t1).^2)).* exp(1i*psy_depth(p1,t1));
    % E_y_muz     = @(p1,t1) -n^2 .*  (cos(t1)./cos2(p1,t1)).*t_p(p1,t1).*rho(p1,t1).*sin(p1).* exp(1i*psy_depth(p1,t1));
    
    E_45_mux    = @(p1,t1) ( E_x_mux(p1,t1) + E_y_mux(p1,t1) )/sqrt(2);
    E_45_muy    = @(p1,t1) ( E_x_muy(p1,t1) + E_y_muy(p1,t1) )/sqrt(2);
    E_45_muz    = @(p1,t1) ( E_x_muz(p1,t1) + E_y_muz(p1,t1) )/sqrt(2);
    
    E_135_mux   = @(p1,t1) (-E_x_mux(p1,t1) + E_y_mux(p1,t1) )/sqrt(2);
    E_135_muy   = @(p1,t1) (-E_x_muy(p1,t1) + E_y_muy(p1,t1) )/sqrt(2);
    E_135_muz   = @(p1,t1) (-E_x_muz(p1,t1) + E_y_muz(p1,t1) )/sqrt(2);
    
    % Factors of the matrix
    
    XXx     = @(p1,t1) abs( E_x_mux(p1,t1) ).^2 .*sin(t1);
    XXy     = @(p1,t1) abs( E_y_mux(p1,t1) ).^2 .*sin(t1);
    
    YYx     = @(p1,t1) abs( E_x_muy(p1,t1) ).^2 .*sin(t1);
    YYy     = @(p1,t1) abs( E_y_muy(p1,t1) ).^2 .*sin(t1);
    
    ZZx     = @(p1,t1) abs( E_x_muz(p1,t1) ).^2 .*sin(t1);
    ZZy     = @(p1,t1) abs( E_y_muz(p1,t1) ).^2 .*sin(t1);
    
    XX45    = @(p1,t1) abs( E_45_mux(p1,t1) ).^2  .*sin(t1);
    XX135   = @(p1,t1) abs( E_135_mux(p1,t1) ).^2 .*sin(t1);
    
    YY45    = @(p1,t1) abs( E_45_muy(p1,t1) ).^2  .*sin(t1);
    YY135   = @(p1,t1) abs( E_135_muy(p1,t1) ).^2 .*sin(t1);
    
    ZZ45    = @(p1,t1) abs( E_45_muz(p1,t1) ).^2  .*sin(t1);
    ZZ135   = @(p1,t1) abs( E_135_muz(p1,t1) ).^2 .*sin(t1);
    
    XYx     = @(p1,t1) 2*real( conj(E_x_mux(p1,t1))   .* E_x_muy(p1,t1) )   .*sin(t1);
    XYy     = @(p1,t1) 2*real( conj(E_y_mux(p1,t1))   .* E_y_muy(p1,t1) )   .*sin(t1);
    
    XY45    = @(p1,t1) 2*real( conj(E_45_mux(p1,t1))  .* E_45_muy(p1,t1) )  .*sin(t1);
    XY135   = @(p1,t1) 2*real( conj(E_135_mux(p1,t1)) .* E_135_muy(p1,t1) ) .*sin(t1);
    
    % XZx = @(p1,t1) 2*real(conj(E_x_mux(p1,t1)) .* E_x_muz(p1,t1)).*sin(t1);
    % XZy = @(p1,t1) 2*real(conj(E_y_mux(p1,t1)) .* E_y_muz(p1,t1)).*sin(t1);
    % 
    % XZ45 = @(p1,t1) 2*real(conj(E_45_mux(p1,t1)) .* E_45_muz(p1,t1)).*sin(t1);
    % XZ135 = @(p1,t1) 2*real(conj(E_135_mux(p1,t1)) .* E_135_muz(p1,t1)).*sin(t1);
    % 
    % YZx = @(p1,t1) 2*real(conj(E_x_muy(p1,t1)) .* E_x_muz(p1,t1)).*sin(t1);
    % YZy = @(p1,t1) 2*real(conj(E_y_muy(p1,t1)) .* E_y_muz(p1,t1)).*sin(t1);
    % 
    % YZ45 = @(p1,t1) 2*real(conj(E_45_muy(p1,t1)) .* E_45_muz(p1,t1)).*sin(t1);
    % YZ135 = @(p1,t1) 2*real(conj(E_135_muy(p1,t1)) .* E_135_muz(p1,t1)).*sin(t1);
    
    %% Intensity integrals
    if red_ch_select == 1
        XX0 = integral2(XXx,0,2*pi,0,theta_red,'method','iterated');
        YY0 = integral2(YYx,0,2*pi,0,theta_red,'method','iterated');
        ZZ0 = integral2(ZZx,0,2*pi,0,theta_red,'method','iterated');
        XY0 = integral2(XYx,0,2*pi,0,theta_red,'method','iterated');
        %XZ0 = integral2(XZx,0,2*pi,0,theta_1red,'method','iterated');
        %YZ0 = integral2(YZx,0,2*pi,0,theta_1red,'method','iterated');
        
        XX90 = integral2(XXy,0,2*pi,0,theta_red,'method','iterated');
        YY90 = integral2(YYy,0,2*pi,0,theta_red,'method','iterated');
        ZZ90 = integral2(ZZy,0,2*pi,0,theta_red,'method','iterated');
        XY90 = integral2(XYy,0,2*pi,0,theta_red,'method','iterated');
        %XZ90 = integral2(XZy,0,2*pi,0,theta_1red,'method','iterated');
        %YZ90 = integral2(YZy,0,2*pi,0,theta_1red,'method','iterated');
        
        XX_45 = integral2(XX45,0,2*pi,0,theta_max,'method','iterated');
        YY_45 = integral2(YY45,0,2*pi,0,theta_max,'method','iterated');
        ZZ_45 = integral2(ZZ45,0,2*pi,0,theta_max,'method','iterated');
        XY_45 = integral2(XY45,0,2*pi,0,theta_max,'method','iterated');
        % XZ_45 = integral2(XZ45,0,2*pi,0,theta_1c,'method','iterated');
        % YZ_45 = integral2(YZ45,0,2*pi,0,theta_1c,'method','iterated');
        
        XX_135 = integral2(XX135,0,2*pi,0,theta_max,'method','iterated');
        YY_135 = integral2(YY135,0,2*pi,0,theta_max,'method','iterated');
        ZZ_135 = integral2(ZZ135,0,2*pi,0,theta_max,'method','iterated');
        XY_135 = integral2(XY135,0,2*pi,0,theta_max,'method','iterated');
        % XZ_135 = integral2(XZ135,0,2*pi,0,theta_1c,'method','iterated');
        % YZ_135 = integral2(YZ135,0,2*pi,0,theta_1c,'method','iterated');
    elseif red_ch_select == 2
            XX0 = integral2(XXx,0,2*pi,0,theta_max,'method','iterated');
            YY0 = integral2(YYx,0,2*pi,0,theta_max,'method','iterated');
            ZZ0 = integral2(ZZx,0,2*pi,0,theta_max,'method','iterated');
            XY0 = integral2(XYx,0,2*pi,0,theta_max,'method','iterated');
            %XZ0 = integral2(XZx,0,2*pi,0,theta_max,'method','iterated');
            %YZ0 = integral2(YZx,0,2*pi,0,theta_max,'method','iterated');
            
            XX90 = integral2(XXy,0,2*pi,0,theta_max,'method','iterated');
            YY90 = integral2(YYy,0,2*pi,0,theta_max,'method','iterated');
            ZZ90 = integral2(ZZy,0,2*pi,0,theta_max,'method','iterated');
            XY90 = integral2(XYy,0,2*pi,0,theta_max,'method','iterated');
            %XZ90 = integral2(XZy,0,2*pi,0,theta_max,'method','iterated');
            %YZ90 = integral2(YZy,0,2*pi,0,theta_max,'method','iterated');
            
            XX_45 = integral2(XX45,0,2*pi,0,theta_red,'method','iterated');
            YY_45 = integral2(YY45,0,2*pi,0,theta_red,'method','iterated');
            ZZ_45 = integral2(ZZ45,0,2*pi,0,theta_red,'method','iterated');
            XY_45 = integral2(XY45,0,2*pi,0,theta_red,'method','iterated');
            % XZ_45 = integral2(XZ45,0,2*pi,0,theta_1c,'method','iterated');
            % YZ_45 = integral2(YZ45,0,2*pi,0,theta_1c,'method','iterated');
            
            XX_135 = integral2(XX135,0,2*pi,0,theta_red,'method','iterated');
            YY_135 = integral2(YY135,0,2*pi,0,theta_red,'method','iterated');
            ZZ_135 = integral2(ZZ135,0,2*pi,0,theta_red,'method','iterated');
            XY_135 = integral2(XY135,0,2*pi,0,theta_red,'method','iterated');
            % XZ_135 = integral2(XZ135,0,2*pi,0,theta_red,'method','iterated');
            % YZ_135 = integral2(YZ135,0,2*pi,0,theta_red,'method','iterated');
    else
    end
    
    K=[XX0 YY0 ZZ0 XY0;...
        XX90 YY90 ZZ90 XY90;...
        XX_45 YY_45 ZZ_45 XY_45;...
        XX_135 YY_135 ZZ_135 XY_135];

end