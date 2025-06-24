function [fit_results, old_data] = fit_4P_ISOmodel(K, varargin)
    %% fit_results = fit_4P_ISOmodel(K, varargin) 
    % Finds dipole orientation parameters:
    %   rho     - in-plane orientation (polar angle)
    %   eta     - out-of-plane orientation (azimuth angle)
    %   delta   - wobbling
    %   
    % Inputs:
    %   K-matrix
    %   variable_intput: 4x4 file from 4polar detection
    % 
    % Output:
    %   fit_results - structure (rho, eta, delta, sqrt(objective function value))
    %       dimension 1: individual detections
    %       dimension 2: "Method"   1-old method without fitting
    %                               2-fit using fmincon initial value real part of eta and delta
    %                               3-fit using fmincon initial value mean of real eta and delta
    % Requires subfunctions: 
    %   - calc_4P_Polar_params_ISO(4x4, K)
    %   - lambda_iso(delta)
    %   - calculate_int(K, rho, eta, delta)
    %   - M_iso(rho, eta, delta)
    %
    % Authors:
    %   Charitra S. Senthil Kumar - Institut Fresnel  
    %   Miguel Sison              - Institut Fresnel
    %
    % Date: March 2024
    %%

    if nargin == 1
        [fname, fpath] = uigetfile();
        file_4x4 = fullfile(fpath,fname);
    elseif nargin == 2
        file_4x4 = varargin{1};
    else
        errordlg('Wrong Inputs!!!')
        return
    end

    data = calc_4P_Polar_params_ISO(file_4x4, K);
    N = numel(data.It);
    
    rho_init    = data.rho;
    eta_init    = data.eta;
    delta_init  = data.delta;
%     x_loc       = data.X(:,1)*0.13;
%     y_loc       = data.Y(:,1)*0.13;
%     time        = data.T;

    I_all       = data.I(:,[1 3 2 4])./data.It; % reorganize 4P intensities columns to 0 90 45 135

    polar_params = NaN(N, 4, 2);
    
    p_init2 = mean([delta_init(~data.ex,1) eta_init(~data.ex,1)], 1);
    temp    = std([delta_init(~data.ex,1) eta_init(~data.ex,1)], [], 1);

    p_max2  = min([p_init2 + 4*temp; pi pi/2], [], 1);
    p_min2  = max([p_init2 - 4*temp; 0 0], [], 1);
    clear temp



    disp(['Fitting orientation paramters... ' num2str(N) ' total detections.'])

    opt = optimoptions(@fmincon);
    opt.Display = 'off';

    delete(gcp('nocreate'));
    % parpool('LocalProfile1');
    parpool('Threads'); % use this if LocalProfile1 doesnt exist or is invalid
    disp('...')
    disp('...')
    disp('...')
    disp('START!!!')
    tic
    parfor k = 1:N
        
        temp_I = I_all(k,:)';
    
        rho_temp    = rho_init(k, 1);
        eta_temp    = real( eta_init(k, 1) );
        delta_temp  = real( delta_init(k, 1) );
    
        %%% Objective function %%%
        fun_temp = @(p) mean( ( temp_I - calculate_int( K, rho_temp, p(2), p(1) ) ).^2, 'all' );
        
        p_init1 = [delta_temp eta_temp];
    
        [p_res1, fval1]   = fmincon( fun_temp,        p_init1, [], [], [], [], [0 0], [pi pi/2], [], opt);
        [p_res2, fval2]   = fmincon( fun_temp,        p_init2, [], [], [], [], p_min2, p_max2, [], opt);
        
%         temp = mean( ( temp_I - calculate_int( K, rho_temp, eta_temp, delta_temp ) ).^2, 'all' );
%         p0  = [rho_temp, eta_temp, delta_temp sqrt( temp )];

        p1  = [rho_temp p_res1([2,1]) sqrt( fval1 )];
        p2  = [rho_temp p_res2([2,1]) sqrt( fval2 )];
        
        polar_params(k,:,:) = cat(3, p1, p2);
    end
    total_sec = toc;
    toc
    delete(gcp('nocreate'));

    temp = sqrt(mean( (I_all' - calculate_int(K, rho_init' , real(eta_init') , real(delta_init'))).^2, 1));
    fit_results.rho     = [rho_init     squeeze( polar_params(:,1,:) )];
    fit_results.eta     = [eta_init     squeeze( polar_params(:,2,:) )];
    fit_results.delta   = [delta_init   squeeze( polar_params(:,3,:) )];
    fit_results.rmse    = [temp'        squeeze( polar_params(:,4,:) )];
    %% XYT section
%     fit_results.X       = [x_loc];
%     fit_results.Y       = [y_loc];
%     fit_results.T       = [time];
%     fit_results.old_data= data;
    %%
    old_data = data;
    disp('DONE!!!')
    disp(['Fitted ' num2str(N) ' total detections in ' num2str(total_sec) ' seconds...'])
end
