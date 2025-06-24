function [Kzz, Kzz_std, I_beads] = Kzz_calib(K_test, K_theo, varargin)
    %% Caculate the Kzz column following the equation
    % I = K*M or Kzz = I - (Kxx*Mxx + Kyy*Myy + Kxy*Mxy)/Mzz
    % and for an ideal isotropic emitter Mxx = Myy = Mzz and Mxy = 0
    %
    % INPUTS
    % K_test    - calibrated K matrix with no Kzz
    % K_theo    - corresponding theoretical K matrix based on bfp reduction
    %             factor
    %
    % varargin  - 4P intensities of isotropic emitter (e.g. fluorescent
    %             beads). If not specified, user will be asked for the 
    %             4x4 file
    %
    % OUTPUT
    % Kzz       - calibrated Kzz column
    %
    % Authors:
    %   Charitra S. Senthil Kumar - Institut Fresnel  
    %   Miguel Sison              - Institut Fresnel
    %
    % Date: March 2024
    %
    %%
    if nargin == 3
        I_beads = varargin{1};
    elseif nargin == 2
        [fname, fpath] = uigetfile();
        beads_file = fullfile(fpath,fname);
        load(beads_file, 'tab_traj_4x4',...
        'param_t_r0', 'param_t_alpha','param_t_tps');
        
        I0      =   (2*sqrt(pi))    .*  tab_traj_4x4(param_t_r0, 1:4:end)    .*  tab_traj_4x4(param_t_alpha, 1:4:end);
        I45     =   (2*sqrt(pi))    .*  tab_traj_4x4(param_t_r0, 2:4:end)    .*  tab_traj_4x4(param_t_alpha, 2:4:end);
        I90     =   (2*sqrt(pi))    .*  tab_traj_4x4(param_t_r0, 3:4:end)    .*  tab_traj_4x4(param_t_alpha, 3:4:end);
        I135    =   (2*sqrt(pi))    .*  tab_traj_4x4(param_t_r0, 4:4:end)    .*  tab_traj_4x4(param_t_alpha, 4:4:end);

        I_beads   = [I0; I90; I45; I135];
        
    else
        errordlg('Wrong inputs')
        return
    end
    % Moment of isotropic emitter
    MM = M_iso(randn,randn,pi);

    % 4P intensities of ideal isotropic emitter and ideal K matrix. because
    % of finite NA the intensities don't normalize to 1
    I_iso = K_theo*MM;
    iso_int_factor = sum(I_iso);

    % normalize isotropic beads intensity and apply the "factor"
    I_beads_norm = I_beads./sum(I_beads)*iso_int_factor;
    
    % resize MM to be the same size as I_beads
    MM = repmat( MM, 1,size(I_beads,2) );
    
    % just to make sure that there is no Kzz term
    % only in plane components
    K_test(:,3) = 0; 
    I_in_plane = K_test*MM;

    % Kzz = I - (Kxx*Mxx + Kyy*Myy + Kxy*Mxy)/Mzz
    % This equation is applied over many fluorescent bead detections and
    % ideally the resulting Kzz should be the same. In this case we take
    % the mean.
    Kzz(:,1)        = mean( (I_beads_norm - I_in_plane)*3 , 2 );
   
    % look at standard deviation
    Kzz_std(:,1)    = std( (I_beads_norm - I_in_plane)*3 , [], 2 );

%     % TRIAL removing outliers
%     clear exc_mat
%     for k = 1:4
%         [~, temp] = rmoutliers(I_beads_norm(k,:), 2);
%         exc_mat(k,:) = temp;
%     end
% 
%     
%     ex = sum(exc_mat,1) < 1;
% 
%     Kzz(:,2)        = mean( (I_beads_norm(:,ex) - I_in_plane(:,ex))*3 , 2 );
%    
%     % look at standard deviation
%     Kzz_std(:,2)    = std( (I_beads_norm(:,ex) - I_in_plane(:,ex))*3 , [], 2 );
   
end