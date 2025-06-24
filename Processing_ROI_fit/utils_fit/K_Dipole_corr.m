function [K_corr, corr_factor] = K_Dipole_corr(K_test, bfp_factor, n_glass, n_sample, NA, SAF_flag, varargin)
    
    if nargin == 6
        lambda      = 0.600;        % wavelength (in um)
        z_dipole    = 0.0; 
    elseif nargin == 7
        lambda      = 0.600;        % wavelength (in um)
        z_dipole    = varargin{1};
    elseif nargin == 8
        lambda      = varargin{2};  % wavelength (in um)
        z_dipole    = varargin{1};
    else
        errordlg('Wrong Inputs!!!')
        return
    end

    % Looks at K_test and checks which channels are reduced. 
    % 0/90 reduced: NAred_ch = 1
    % 45/135 reduced: NAred_ch = 2
    [~, NAred_ch] = min( mean(K_test(1:2:end,1:2) + K_test(2:2:end,1:2),2) );
       
    %% set max NA - with (full objective NA) or without SAF (water glass interface)
    
    n = n_glass/n_sample;
    
    if SAF_flag
        theta_max = asin(NA./n_glass);
        theta_red = asin(bfp_factor.*NA./n_glass);
    else
        theta_max = asin(n_sample./n_glass);
        theta_red = asin(bfp_factor.*n_sample./n_glass);
    end

    %% General parameters and functions for integration
    
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
    XX000   = @(p1,t1) abs(E_x_mux(p1,t1)).^2 .* sin(t1);
    XX090   = @(p1,t1) abs(E_y_mux(p1,t1)).^2 .* sin(t1);
    XX045   = @(p1,t1) abs(E_45_mux(p1,t1)).^2 .* sin(t1);
    XX135   = @(p1,t1) abs(E_135_mux(p1,t1)).^2 .* sin(t1);

    YY000   = @(p1,t1) abs(E_x_muy(p1,t1)).^2 .* sin(t1);
    YY090   = @(p1,t1) abs(E_y_muy(p1,t1)).^2 .* sin(t1);
    YY045   = @(p1,t1) abs(E_45_muy(p1,t1)).^2 .* sin(t1);
    YY135   = @(p1,t1) abs(E_135_muy(p1,t1)).^2 .* sin(t1);

    ZZ000   = @(p1,t1) abs(E_x_muz(p1,t1)).^2  .* sin(t1);
    ZZ090   = @(p1,t1) abs(E_y_muz(p1,t1)).^2  .* sin(t1);   
    ZZ045   = @(p1,t1) abs(E_45_muz(p1,t1)).^2 .* sin(t1);
    ZZ135   = @(p1,t1) abs(E_135_muz(p1,t1)).^2 .* sin(t1);
    
    XY000   = @(p1,t1) 2*real(conj(E_x_mux(p1,t1)) .* E_x_muy(p1,t1)) .* sin(t1);
    XY090   = @(p1,t1) 2*real(conj(E_y_mux(p1,t1)) .* E_y_muy(p1,t1)) .* sin(t1);
    XY045   = @(p1,t1) 2*real(conj(E_45_mux(p1,t1)) .* E_45_muy(p1,t1)) .* sin(t1);
    XY135   = @(p1,t1) 2*real(conj(E_135_mux(p1,t1)) .* E_135_muy(p1,t1)) .* sin(t1);
    
    XZ000   = @(p1,t1) 2*real(conj(E_x_mux(p1,t1)) .* E_x_muz(p1,t1)) .* sin(t1);
    XZ090   = @(p1,t1) 2*real(conj(E_y_mux(p1,t1)) .* E_y_muz(p1,t1)) .* sin(t1);
    XZ045   = @(p1,t1) 2*real(conj(E_45_mux(p1,t1)) .* E_45_muz(p1,t1)) .* sin(t1);
    XZ135   = @(p1,t1) 2*real(conj(E_135_mux(p1,t1)) .* E_135_muz(p1,t1)) .* sin(t1);
    
    YZ000   = @(p1,t1) 2*real(conj(E_x_muy(p1,t1)) .* E_x_muz(p1,t1)) .* sin(t1);
    YZ090   = @(p1,t1) 2*real(conj(E_y_muy(p1,t1)) .* E_y_muz(p1,t1)) .* sin(t1);
    YZ045   = @(p1,t1) 2*real(conj(E_45_muy(p1,t1)) .* E_45_muz(p1,t1)) .* sin(t1);
    YZ135   = @(p1,t1) 2*real(conj(E_135_muy(p1,t1)) .* E_135_muz(p1,t1)) .* sin(t1);
      
    %% theoretical K matrix from ideal rotating dipole

    if NAred_ch == 1
        % intensity integrals high NA channels
       
        Kxx045 = integral2(XX045,0,2*pi,0,theta_max,'method','iterated');
        Kyy045 = integral2(YY045,0,2*pi,0,theta_max,'method','iterated');
        Kzz045 = integral2(ZZ045,0,2*pi,0,theta_max,'method','iterated');
        Kxy045 = integral2(XY045,0,2*pi,0,theta_max,'method','iterated');
        Kxz045 = integral2(XZ045,0,2*pi,0,theta_max,'method','iterated');
        Kyz045 = integral2(YZ045,0,2*pi,0,theta_max,'method','iterated');
        
        Kxx135 = integral2(XX135,0,2*pi,0,theta_max,'method','iterated');
        Kyy135 = integral2(YY135,0,2*pi,0,theta_max,'method','iterated');
        Kzz135 = integral2(ZZ135,0,2*pi,0,theta_max,'method','iterated');
        Kxy135 = integral2(XY135,0,2*pi,0,theta_max,'method','iterated');
        Kxz135 = integral2(XZ135,0,2*pi,0,theta_max,'method','iterated');
        Kyz135 = integral2(YZ135,0,2*pi,0,theta_max,'method','iterated');
    
        % intensity integrals low NA channels
        Kxx000 = integral2(XX000,0,2*pi,0,theta_red,'method','iterated');
        Kyy000 = integral2(YY000,0,2*pi,0,theta_red,'method','iterated');
        Kzz000 = integral2(ZZ000,0,2*pi,0,theta_red,'method','iterated');
        Kxy000 = integral2(XY000,0,2*pi,0,theta_red,'method','iterated');
        Kxz000 = integral2(XZ000,0,2*pi,0,theta_red,'method','iterated');
        Kyz000 = integral2(YZ000,0,2*pi,0,theta_red,'method','iterated');
        
        Kxx090 = integral2(XX090,0,2*pi,0,theta_red,'method','iterated');
        Kyy090 = integral2(YY090,0,2*pi,0,theta_red,'method','iterated');
        Kzz090 = integral2(ZZ090,0,2*pi,0,theta_red,'method','iterated');
        Kxy090 = integral2(XY090,0,2*pi,0,theta_red,'method','iterated');
        Kxz090 = integral2(XZ090,0,2*pi,0,theta_red,'method','iterated');
        Kyz090 = integral2(YZ090,0,2*pi,0,theta_red,'method','iterated');

    elseif NAred_ch == 2
        % intensity integrals high NA channels
        Kxx000 = integral2(XX000,0,2*pi,0,theta_max,'method','iterated');
        Kyy000 = integral2(YY000,0,2*pi,0,theta_max,'method','iterated');
        Kzz000 = integral2(ZZ000,0,2*pi,0,theta_max,'method','iterated');
        Kxy000 = integral2(XY000,0,2*pi,0,theta_max,'method','iterated');
        Kxz000 = integral2(XZ000,0,2*pi,0,theta_max,'method','iterated');
        Kyz000 = integral2(YZ000,0,2*pi,0,theta_max,'method','iterated');
        
        Kxx090 = integral2(XX090,0,2*pi,0,theta_max,'method','iterated');
        Kyy090 = integral2(YY090,0,2*pi,0,theta_max,'method','iterated');
        Kzz090 = integral2(ZZ090,0,2*pi,0,theta_max,'method','iterated');
        Kxy090 = integral2(XY090,0,2*pi,0,theta_max,'method','iterated');
        Kxz090 = integral2(XZ090,0,2*pi,0,theta_max,'method','iterated');
        Kyz090 = integral2(YZ090,0,2*pi,0,theta_max,'method','iterated');

        % intensity integrals low NA channels
        Kxx045 = integral2(XX045,0,2*pi,0,theta_red,'method','iterated');
        Kyy045 = integral2(YY045,0,2*pi,0,theta_red,'method','iterated');
        Kzz045 = integral2(ZZ045,0,2*pi,0,theta_red,'method','iterated');
        Kxy045 = integral2(XY045,0,2*pi,0,theta_red,'method','iterated');
        Kxz045 = integral2(XZ045,0,2*pi,0,theta_red,'method','iterated');
        Kyz045 = integral2(YZ045,0,2*pi,0,theta_red,'method','iterated');
        
        Kxx135 = integral2(XX135,0,2*pi,0,theta_red,'method','iterated');
        Kyy135 = integral2(YY135,0,2*pi,0,theta_red,'method','iterated');
        Kzz135 = integral2(ZZ135,0,2*pi,0,theta_red,'method','iterated');
        Kxy135 = integral2(XY135,0,2*pi,0,theta_red,'method','iterated');
        Kxz135 = integral2(XZ135,0,2*pi,0,theta_red,'method','iterated');
        Kyz135 = integral2(YZ135,0,2*pi,0,theta_red,'method','iterated');
    else
        errordlg('Wrong choice for NAred_ch! Choose either 1 (0/90 reduced NA) or 2 (45/135 reduced NA)');
        return
    end
    
    K_theo       = [Kxx000 Kyy000 Kzz000 Kxy000; 
                   Kxx090 Kyy090 Kzz090 Kxy090; 
                   Kxx045 Kyy045 Kzz045 Kxy045; 
                   Kxx135 Kyy135 Kzz135 Kxy135];
    K_theo       = K_theo./sum(K_theo(:,1));

%     K_theo_full  = [Kxx000 Kyy000 Kzz000 Kxy000 Kxz000 Kyz000; 
%                    Kxx090 Kyy090 Kzz090 Kxy090 Kxz090 Kyz090; 
%                    Kxx045 Kyy045 Kzz045 Kxy045 Kxz045 Kyz045; 
%                    Kxx135 Kyy135 Kzz135 Kxy135 Kxz135 Kyz135];
%     
%     K_theo_full  = K_theo_full./sum(K_theo_full(:,1));


    %% theoretical intensities on all 4-polar channels
    
    N_steps = 180;
    alpha = deg2rad(linspace(0,180,N_steps+1));
    alpha = alpha(1:end-1);
    
    Mu      = [cos(alpha).^2;
               sin(alpha).^2;
               zeros(size(alpha));
               cos(alpha).*sin(alpha)];

%     Mu_full = [cos(alpha).^2;
%                sin(alpha).^2;
%                zeros(size(alpha));
%                cos(alpha).*sin(alpha);
%                zeros(size(alpha)).*cos(alpha);
%                zeros(size(alpha)).*sin(alpha)];
%     
%     I_mu_full    = K_theo_full*Mu_full; 
%     I_mu = I_mu_full;

    % rotating dipole
    I_mu    = K_theo*Mu; % already normalized only if K_theo and Mu are 
               
    % rotating polarizer
    I_pol   = sum(K_theo(:,1:3),2) .* cos(alpha+[0; pi/2; -pi/4; +pi/4]).^2 ;
    % I_pol   = [Kxx000 + Kyy000 + Kzz000; Kxx090 + Kyy090 + Kzz090; Kxx045 + Kyy045 + Kzz045; Kxx135 + Kyy135 + Kzz135] .* cos(alpha+[0; pi/2; -pi/4; +pi/4]).^2 ;
    % these 2 expressions for I_pol are the same

%     % normalize
%     I_pol   = I_pol./sum(I_pol, 1);
%     I_mu    = I_mu./sum(I_mu, 1);
    
    % Calculate even Fourier series components (up to n = 2) of I_mu and I_pol (circular decomposition)
    I_temp  = I_mu;
    I_temp  = I_temp./sum(I_temp,1);
    A0mu    = mean( I_temp, 2 );
    A2mu    = 2*mean( I_temp.*cos(2*alpha), 2 );
    B2mu    = 2*mean( I_temp.*sin(2*alpha), 2 );
       
    I_temp  = I_pol;
    I_temp  = I_temp./sum(I_temp,1);
    A0pol   = mean( I_temp, 2 );
    A2pol   = 2*mean( I_temp.*cos(2*alpha), 2 );
    B2pol   = 2*mean( I_temp.*sin(2*alpha), 2 );

    clear I_temp

    % Calculate correction factors with normalization. 
    % zero_adj is for a correcting or avoiding dividing by zero. A better 
    % solution might be needed, taking limit for example
    zero_adj = 0.00000001;
    zero_adj = max(100000*min(abs([A0mu A2mu B2mu A0pol A2pol B2pol]),[],'all'), zero_adj);
%     zero_adj = 0.00000;

%     corr_factor = ( [A0mu A2mu B2mu] + zero_adj.*sign([A0mu A2mu B2mu]) )...
%                     ./   ...
%                   ( [A0pol A2pol B2pol] + zero_adj.*sign([A0pol A2pol B2pol]) );
    corr_factor = ( max(abs([A0mu A2mu B2mu]),zero_adj) )...
                    ./   ...
                  ( max(abs([A0pol A2pol B2pol]),zero_adj) );

    corr_factor = abs(corr_factor);
    %% Generate "experimental intensities" based on K_test and calculate 
    % even Fourier series components up to n = 2

%     % Use if you want different number of test angles
%     N = 180;
%     test_alpha = linspace(0,pi,N+1);
%     test_alpha = test_alpha(1:end-1);
%     
%     Mu_test    = [cos(test_alpha).^2;
%                   sin(test_alpha).^2;
%                   zeros(size(test_alpha));
%                   cos(test_alpha).*sin(test_alpha)];
%     I_test      = K_all*Mu_test;

    test_alpha = alpha;

    I_test  = K_test*Mu;
    I_test  = I_test./sum(I_test,1);
    
    I_temp  = I_test;
    
    A0test  = mean( I_temp, 2 );
    A2test  = 2*mean( I_temp.*cos(2*test_alpha), 2 );
    B2test  = 2*mean( I_temp.*sin(2*test_alpha), 2 );

    % Apply correction factors to components
    temp = [A0test A2test B2test] .* corr_factor;
    A0corr = temp(:,1);
    A2corr = temp(:,2);
    B2corr = temp(:,3);

%     [A0mu A2mu B2mu]
%     [A0pol A2pol B2pol]
%     [A0test A2test B2test]
%     [A0corr A2corr B2corr]
    
    % generate new intensities A0 + A2*cos(2*alpha) + B2*sin(2*alpha)
    I_corr  = temp*[ones(size(test_alpha));
                    cos(2*test_alpha);
                    sin(2*test_alpha)];
    clear temp
    %% Calculate corrected K
    
%     [inds,~] = find((test_alpha == [0; pi/2; pi/4; 3*pi/4])');
% 
%     temp = I_corr(:, inds);
%     K_corr(:,1:2)    = temp(:,1:2);
%     K_corr(:,4)      = 2*temp(:,3) - sum(temp(:,1:2),2); %2*I(45) - I(0) - I(90)
%     K_corr = K_corr./sum(K_corr(:,1));
    
    % % % I = A0 + A2*cos(2*alpha) + B2*sin(2*alpha)
    % % % I = (A0 + A2)*cos(alpha).^2 + (A0 - A2)*sin(alpha).^2 + 2*B2*sin(alpha)*cos(alpha)
    K_corr(:, 1:2)    = [A0corr + A2corr, A0corr - A2corr];
    K_corr(:,4)       = 2*B2corr;
    K_corr            = K_corr./sum(K_corr(:,1),1);
   
    K_corr(:,3)      = K_test(:,3); % Kzz from K_test
    %% look at K's
    % K_mu:     from rotating dipole
    % K_pol:    from rotating polarizer
    
    % [inds,~] = find((alpha == [0;pi/2;pi/4;3*pi/4])');
    % 
    % temp = I_mu(:,inds);
    % 
    % K_mu(:,1:2)     = temp(:,1:2);
    % K_mu(:,4)       = 2*temp(:,3) - sum(temp(:,1:2),2); %2*I(45) - I(0) - I(90)
    % 
    % K_mu(:,3)       = K_theo(:,3); % theoretical Kzz
    % K_mu = K_mu./sum(K_mu(:,1));
    % 
    % temp = I_pol(:,inds);
    % 
    % K_pol(:,1:2)    = temp(:,1:2);
    % K_pol(:,4)      = 2*temp(:,3) - sum(temp(:,1:2),2); %2*I(45) - I(0) - I(90)
    % 
    % K_pol(:,3)       = K_theo(:,3); % theoretical Kzz
    % K_pol = K_pol./sum(K_pol(:,1));
    %
    % or 
    %
    % K_mu(:, 1:2)    = [A0mu + A2mu, A0mu - A2mu];
    % K_mu(:,4)       = 2*B2mu;
    % K_mu            = K_mu./sum(K_mu(:,1),1);
    % 
    % K_pol(:, 1:2)    = [A0pol + A2pol, A0pol - A2pol];
    % K_pol(:,4)       = 2*B2pol;
    % K_pol            = K_pol./sum(K_pol(:,1),1);

    
end