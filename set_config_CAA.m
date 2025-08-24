function cfg = set_config_CAA(config_tag)
% Configuration settings for the CAA dataset
clear cfg;

cfg.only_post_processing    = 1; % set to 1 if already run ROMT, 
% this could skip ROMT and only do GLAD using serialized u and rho fields.
cfg.dataset_name            = 'CAA';  % name of dataset (do not change, fix as standard setting)

cfg.tag                     = config_tag; % name result folder (could change to any name)

config_tag_prefix = split(config_tag, '_');
config_tag_postfix = config_tag_prefix{end}; % get the last part of the tag, e.g., 'test1' or 'C294'
config_tag_prefix = config_tag_prefix{1};

switch config_tag_prefix

case 'ours'
    % preprocessed matdata, if set, will cover nii data.
    % cfg.matdata_template        = './data/ours/test1/DCE_mat_data/dce%d.mat'; % template for mat data loading,
    % cfg.matdata_fieldname       = 'V';
    switch config_tag_postfix
        case 'test1'
            cfg.data_template           = './data/ours/test1/DCE_nii_data/psnrT1_FLASH_3D_%04d.nii'; % template for data loading
            cfg.ROI_msk_path            = './data/ours/test1/Template_C57Bl6_n30_brain.nii'; % basic mask to set control volume.
        case 'DEXI' 
            cfg.data_template           = '/data/xym/DEX_MRI/DEXI/DEXI_083/DCE_nii_data/psnrT1_FLASH_3D_%04d.nii'; % template for data loading
            cfg.ROI_msk_path            = '/data/xym/DEX_MRI/DEXI/Template_C57Bl6_n30_brain_DEXI_083.nii'; % basic mask to set control volume.
            cfg.exclude_frames          = [4, 5, 6];
        case 'ISO'
            cfg.data_template           = '/data/xym/DEX_MRI/ISO/ISO_52/DCE_nii_data/psnrT1_FLASH_3D_%04d.nii'; % template for data loading
            cfg.ROI_msk_path            = '/data/xym/DEX_MRI/ISO/Template_C57Bl6_n30_brain_ISO_52.nii'; % basic mask to set control volume.
            cfg.exclude_frames          = [4, 5, 6, 7, 8];
            
            % WARNING: A Nx3x3 tensor field, N=XYZ. If DTI path is set, then diffuse term (laplace matrix) will calculate using getAnisotropicDiffusion.m
            % and dti_enhanced will be used, sigma will be ignored.
            cfg.dti_path                = '/data/xym/DTI_data/ISO 052/dicom/out_dwi/dti_aligned/dti_tensor.mat';
        case 'KX'
            cfg.data_template           = '/data/xym/DEX_MRI/KX/KX_078/DCE_nii_data/psnrT1_FLASH_3D_%04d.nii'; % template for data loading
            cfg.ROI_msk_path            = '/data/xym/DEX_MRI/KX/Template_C57Bl6_n30_brain_KX_078.nii'; % basic mask to set control volume.
            cfg.exclude_frames          = [4, 5, 6, 7, 8, 27, 29, 30, 31]; % set the speeed field of exclude_frames to zero.

            cfg.dti_path                = '/data/xym/DTI_data/KX 078/dicom/out_dwi/dti_aligned/dti_tensor.mat';
        otherwise
            error('Unknown dataset postfix: %s', config_tag_postfix);
    end

    % basic mask for setting control volume + background anatomy image.
    cfg.do_ROI_msk              = 1;
    cfg.ROI_msk_threshold       = 0.01; % threshold for mask, if msk = raw_msk>0 cause problem.
    
    cfg.x_range                 = 1:128;
    cfg.y_range                 = 1:160;
    cfg.z_range                 = 1:128;

    % set rOMT parameters
    cfg.do_resize               = 1;%1;
    cfg.size_factor             = 0.5;%0.5;
    %cfg.data_index_E            = 7:30;
    
    cfg.smooth                  = 1.2; % evolution time when doing diffusion process runs.Larger t_tot values result in more smoothing. 
    % Think of it as the "exposure time" for the blurring effect.
    
    cfg.dilate                  = 3; % dilate monitor zone, do not miss any volume.
    
    % set bigger if you want to use more data, e.g., 7:50 for 44 time points.
    cfg.first_time              = 7; %9;%cfg.data_index_E(13);，do not include inject time.
    cfg.time_jump               = 1; %3;
    cfg.last_time               = 31;%cfg.data_index_E(33);%;cfg.data_index_E(31);
    
    % empirically set parameters
    cfg.dti_enhanced            = 3.0; % when use dti_path, better enhance the tensor field if Mean Diffusivity too small. 
    cfg.sigma                   = 2e-3; % diffusion coefficient, when cfg.dti_path is set, it will become useless
    cfg.dt                      = 1.0; % 0.2;% timestep for every steps among nt*(last_time-first_time)/time_jump ,
    
    % for DCE-MRI seq, make sure dt * nt = 4 , so that velocity's unit is grid/min
    % or make sure dt * nt = 4/60/0.125 = 0.5333, so that velocity's unit is mm/s
    
    cfg.nt                      = 4; % interpolate velocity field, smaller timestep for easier convergence + memory usage
    
    % WARING: loss of R_1 and R_3 loss is lot smaller than R_2 (phi_match loss),
    cfg.gamma                   = 0.008; % rOMT R_3 loss, for spatial smoothness of velocity field.
    cfg.beta                    = 0.0001; % rOMT R_1 loss, kinetic energy of velocity field, set bigger if not stable.
    
    % if do not have frame_number times memory, replace "parfor" in runROMT_par.m with "par", and also cancel this to get smoother result.
    cfg.reinitR                 = 1; % (No parallel version only) if do consecutively and 1 if reinitialize rho
    cfg.reInitializeU           = 1; % (No parallel version only), 1 if reinitialize u to 0 before each time step; 0 if not, unless first time step
    
    cfg.niter_pcg               = 500; % rounds for pcg solver, set higher if Hl=-g hard to converge. together with update steps par.maxUiter;
    cfg.maxUiter                = 40; % step for update of velocity field, set higher if convergence(see from loss figure) is hard.
    
    cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
    cfg.add_source              = 0; % for unbalanced rOMT.
    
    % GLAD2 config
    % filter for GLAD pathline source points area.
    cfg.density_percent_thres   = 2;  % Threshold to mask velocity, (e.g., 0.1 = 10% bigger than baseline signal)
    cfg.sp_thresh               = 2;  % MRI intensity higher than sp_thresh percentage will consider as start point.
    cfg.GLAD_spfs               = 18;  % starting area sample interval
    cfg.sl_tol                  = 1.2; %threshold for minimum Euclidean length between initial and final streamline points

    % masks for GLAD pathline source points area.
    %cfg.max_dpsnrv             = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii'; %
    cfg.sp_mask_opts(1).name    = 'brain'; % name of the mask, activate it by delete "_disabled". 
    cfg.sp_mask_opts(1).path    =  cfg.ROI_msk_path; %'./data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';
    cfg.sp_mask_opts(1).threshold = 0.01; % for tissue probability mask, 0.5 is standard for CFS/W/G area.
    % update step

    cfg.integral_euler_step     = 10; % number of Eulerian steps to be taken per dt, higher for accurate.
    % trick of extend pathline: fake timestep.
    cfg.GLAD_timestep_factor    = 1; % factor to extend the pathline length, 2.5 means 150% longer pathline.

    % visualization settings
    cfg.speedmap_slice          = 12;
    cfg.anato                   = './data/ours/test1/Template_C57Bl6_n30_brain.nii'; % for vtk viewing, no use.
    cfg.view_azi_elevation      = [5, 80]; % view angle for visualization

    cfg.strid                   = 1; % stride for visualization of flux vectors   
    cfg.flip_z                  = 0; % flip z-axis for visualization
    
    cfg.vis_font_size = 16; % Unified font size for all visualization titles and axes


case 'C294'
    % set directories
    cfg.data_template = './data/12_MONTH_DATA/psnrv/WT/C294/C294_031318A/psnrv_C294_031318A_DOTA37_30ul_E%02d.nii'; 

    cfg.do_ROI_msk              = 1;
    cfg.ROI_msk_path            = './data/12_MONTH_DATA/masksforOMT/WT/C294_MASK.nii';
    cfg.ROI_msk_threshold       = 0.01; % threshold for mask, if msk = raw_msk>0 cause problem.
    cfg.x_range                 = 20:80;
    cfg.y_range                 = 1:106;
    cfg.z_range                 = 39:85;

    %cfg.max_dpsnrv              = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii'; %
    cfg.anato                   = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii'; %
    cfg.sp_mask_opts(1).name    = 'brain';
    cfg.sp_mask_opts(1).path    = './data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';
    cfg.sp_mask_opts(1).threshold = 0; % for tissue probability mask, 0.5 is standard for CFS/W/G area.
    % set rOMT parameters
    cfg.do_resize               = 0;%1;
    cfg.size_factor             = 1;%0.5;
    %cfg.data_index_E            = 19:53;
    cfg.smooth                  = 1;

    cfg.reinitR                 = 1;%1;%0; %0 if do consecutively (runROMT.m) ;
    % and 1 if reinitialize rho using ground truth, and runROMT_par

    cfg.dilate                  = 3;

    cfg.first_time              = 31;%cfg.data_index_E(13);

    cfg.time_jump               = 2;%3;

    cfg.last_time               = 51;%cfg.data_index_E(33);%;cfg.data_index_E(31);

    cfg.sigma                   = 2e-3;
    cfg.dt                      = 0.4; % time step, only affect accuracy because optimization will adjust velocity according to dt

    cfg.nt                      = 4;

    cfg.gamma                   = 0.008;
    cfg.beta                    = 0.0001;
    cfg.reInitializeU           = 1; % (No parallel version only), leverage phi field from last step for better convergence.
    cfg.niter_pcg               = 20;%20;
    cfg.maxUiter                = 25; % step for update of velocity field, set higher if convergence is hard.

    cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
    cfg.add_source              = 0;

    % GLAD2 config
    cfg.GLAD_spfs               = 2; % starting area sample interval.
    cfg.sp_thresh               = 0.2; % higher than sp_thresh percentage will consider as start point.
    cfg.integral_euler_step     = 10; % number of Eulerian steps to be taken per dt, higher for accurate.
    cfg.sl_tol                  = 1; %threshold for minimum Euclidean length between initial and final streamline points
    cfg.GLAD_timestep_factor    = 2.5; % factor to extend the pathline length, 1.5 means 50% longer pathline.
    cfg.density_percent_thres   = 1.05;  % Threshold to mask velocity, (e.g., 0.1 = 10% of baseline signal)

    % visualization settings
    cfg.speedmap_slice          = 10;
    cfg.view_azi_elevation      = [-188.3500   13.7463]; % view angle for visualization
    cfg.strid                   = 10; % stride for visualization of flux vectors   
    cfg.flip_z                  = 0;

    cfg.vis_font_size = 16; % Unified font size for all visualization titles and axes

end

end
