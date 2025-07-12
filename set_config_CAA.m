function cfg = set_config_CAA()
% Configuration settings for the CAA dataset
clear cfg;

cfg.only_post_processing    = 1; % set to 1 if already run ROMT, 
% this could skip ROMT and only do GLAD using serialized u and rho fields.
cfg.dataset_name            = 'CAA';  % name of dataset (do not change, fix as standard setting)

if exist('config_label', 'var') && strcmp(config_label, "ori")
    cfg.tag                     = 'test1'; % name result folder (could change to any name)

    % set directories

    % preprocessed matdata, if set, will cover nii data.
    % cfg.matdata_template        = './data/ours/test1/DCE_mat_data/dce%d.mat'; % template for mat data loading,
    % cfg.matdata_fieldname       = 'V';
    cfg.data_template           = './data/ours/test1/DCE_nii_data/T1_FLASH_3D_%04d.nii'; % template for data loading

    % basic mask for setting control volume + background anatomy image.
    cfg.do_ROI_msk              = 1;
    cfg.ROI_msk_path            = './data/ours/test1/test1_ROImask_new.nii'; % basic mask to set control volume.
    cfg.x_range                 = 1:128;
    cfg.y_range                 = 1:160;
    cfg.z_range                 = 1:128;


    % set rOMT parameters
    cfg.do_resize               = 0;%1;
    cfg.size_factor             = 1;%0.5;
    %cfg.data_index_E            = 7:30;
    cfg.smooth                  = 1; % 总平滑时间（虚拟时间单位，非实际秒数） / dt = 10
    cfg.reinitR                 = 1;%1;%0; %0 if do consecutively and 1 if reinitialize rho
    cfg.dilate                  = 3; % dilate monitor zone, do not miss any volume.

    cfg.first_time              = 4; %9;%cfg.data_index_E(13);，do not include inject time.
    cfg.time_jump               = 1; %3;
    cfg.last_time               = 30;%cfg.data_index_E(33);%;cfg.data_index_E(31);

    % empirically set parameters
    cfg.sigma                   = 2e-3; % diffusion coefficient

    cfg.dt                      = 0.2; % 0.2;% timestep for every steps among nt*(last_time-first_time)/time_jump ,
    % make velocity sensitive by setting dt small.

    cfg.nt                      = 2; % interpolate velocity field, smaller timestep for easier convergence + memory usage

    cfg.gamma                   = 0.008; % rOMT R_3 loss, for smoothness of velocity field.
    cfg.beta                    = 0.001; % rOMT R_1 loss, kinetic energy of velocity field.


    cfg.reInitializeU           = 1; % for better convergence, could leverage speed field in last timestep.
    cfg.niter_pcg               = 30; % rounds for pcg solver, set higher if Hl=-g hard to converge. together with update steps par.maxUiter;
    cfg.maxUiter                = 50; % step for update of velocity field, set higher if convergence is hard.

    cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
    cfg.add_source              = 0; % for unbalanced rOMT.


    % GLAD2 config
    % filter for GLAD pathline source points area.
    cfg.sp_thresh               = 0.05; % higher than sp_thresh percentage will consider as start point.
    cfg.GLAD_spfs               = 2; % starting area sample interval
    cfg.sl_tol                  = 1.5; %threshold for minimum Euclidean length between initial and final streamline points

    % masks for GLAD pathline source points area.
    %cfg.max_dpsnrv              = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii'; %
    cfg.sp_mask_opts(1).name    = 'brain_disabled'; % name of the mask, activate it by delete "_disabled". 
    cfg.sp_mask_opts(1).path    =  './data/ours/test1/test1_spmask_new.nii'; %'./data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';

    % update step
    cfg.integral_euler_step     = 4; % number of Eulerian steps to be taken per dt, higher for accurate.
    % trick of extend pathline: fake timestep.
    cfg.GLAD_timestep_factor    = 2.5; % factor to extend the pathline length, 1.5 means 50% longer pathline.

    % visualization settings
    cfg.speedmap_slice          = 10;
    cfg.anato                   = './data/ours/test1/anatomy_new.nii'; % for vtk viewing, no use.
    cfg.view_azi_elevation      = [20, 70]; % view angle for visualization

    cfg.strid                   = 10; % stride for visualization of flux vectors   
    cfg.flip_z                  = 0; % flip z-axis for visualization

else
    cfg.tag                     = 'C294'; % name of this case
    % set directories
    cfg.data_template = './data/12_MONTH_DATA/psnrv/WT/C294/C294_031318A/psnrv_C294_031318A_DOTA37_30ul_E%02d.nii'; 

    cfg.do_ROI_msk              = 1;
    cfg.ROI_msk_path            = './data/12_MONTH_DATA/masksforOMT/WT/C294_MASK.nii';
    cfg.x_range                 = 20:80;
    cfg.y_range                 = 1:106;
    cfg.z_range                 = 39:85;

    %cfg.max_dpsnrv              = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii'; %
    cfg.anato                   = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii'; %
    cfg.sp_mask_opts(1).name    = 'brain';
    cfg.sp_mask_opts(1).path    = './data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';

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
    cfg.reInitializeU           = 1; % for better convergence, could leverage speed field in last timestep.
    cfg.niter_pcg               = 30;%20;
    cfg.maxUiter                = 50; % step for update of velocity field, set higher if convergence is hard.

    cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
    cfg.add_source              = 0;

    % GLAD2 config
    cfg.GLAD_spfs               = 2; % starting area sample interval.
    cfg.sp_thresh               = 0.1; % higher than sp_thresh percentage will consider as start point.
    cfg.integral_euler_step     = 2; % number of Eulerian steps to be taken per dt, higher for accurate.
    cfg.sl_tol                  = 2; %threshold for minimum Euclidean length between initial and final streamline points
    cfg.GLAD_timestep_factor    = 1.5; % factor to extend the pathline length, 1.5 means 50% longer pathline.

    % visualization settings
    cfg.speedmap_slice          = 10;
    cfg.view_azi_elevation      = [-188.3500   13.7463]; % view angle for visualization
    cfg.strid                   = 10; % stride for visualization of flux vectors   
    cfg.flip_z                  = 0;
end

end
