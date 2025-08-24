% This functions contains the rOMT algorithm and post-processing of the
% results on a sample dataset named 'C294' which is from a CAA dataset containsing about 48 rat brain cases w/w.o CAA disease.

addpath('./Inverse','./Sensitivities','./analyzeFlows',genpath('./utilities'))

%% set directories and parameters

if ~exist('config_tag', 'var')
    config_tag = 'ours_test1'; % default label, could be changed to "C294" for C294 dataset.
end

cfg = set_config_CAA(config_tag); % Load configuration from function

if ~isfield(cfg, 'exclude_frames')
    cfg.exclude_frames = [];
end

if isfield(cfg, 'dti_path') && isfield(cfg, 'stagger_dti_path')
    cfg.D_tensor = load(cfg.dti_path, 'D_tensor').D_tensor;
    if cfg.do_resize
        cfg.D_tensor = resizeDTIMatrix(cfg.D_tensor, round(cfg.size_factor.*size(cfg.D_tensor,1:3)),'linear');
    end
    % cfg.stagger_D_tensor = load(cfg.stagger_dti_path, 'stagger_D_tensor').stagger_D_tensor;
    fprintf("Loaded DTI tensor from: %s\n", cfg.dti_path);
else
    cfg.D_tensor = []; % to trigger isempty and get isotropic laplace matrix in tryGetAnisotropicDiffusion.m
end

% load ROI
if cfg.do_ROI_msk
    tmp = nii2mat(cfg.ROI_msk_path,cfg.x_range,cfg.y_range,cfg.z_range);
    cfg.msk = tmp>1e-4; % WARNING: adjust threshold if necessary
else
    cfg.msk = ones(length(cfg.x_range),length(cfg.y_range),length(cfg.z_range));
end
if cfg.do_resize
   cfg.msk = resizeMatrix(cfg.msk,round(cfg.size_factor.*size(cfg.msk)),'linear');
   cfg.msk(cfg.msk~=1) = 0;
end
if cfg.dilate>0
    [xr,yr,zr] = meshgrid(-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate);
    strel = (xr/cfg.dilate).^2 + (yr/cfg.dilate).^2 + (zr/cfg.dilate).^2 <= 1;
    cfg.msk = imdilate(cfg.msk,strel);
end



% load vol

global_steps = (cfg.last_time-cfg.first_time)/cfg.time_jump+2;
vol_tmp = cell(global_steps, 1);

num_iter = global_steps;

% if has parpool, not allocate again
if ~isempty(gcp('nocreate'))
    num_iter = gcp().NumWorkers; % Get the number of workers in the current parallel pool
else
    parpool(num_iter); % Or parpool(12), etc.
end

parfor i = 1:global_steps
    cur_frame = cfg.first_time+(i-1)*cfg.time_jump;
    fprintf('Loading frame %d', cur_frame);
    % WARNING: change from 02d to 04d if for different dataloader.
    if isfield(cfg, "matdata_template")
        tmp = load(sprintf(cfg.matdata_template, cur_frame));
        tmp = tmp.(cfg.matdata_fieldname);
    else
        tmp = nii2mat(sprintf(cfg.data_template, cur_frame),cfg.x_range,cfg.y_range,cfg.z_range);
    end
    if cfg.do_resize
       tmp = resizeMatrix(tmp,round(cfg.size_factor.*size(tmp)),'linear');
    end
    if cfg.smooth>0
        tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    tmp(~cfg.msk) = 0;
    vol_tmp{i} = tmp;
end


cfg.sig_str                 = erase(num2str(cfg.sigma,'%.0e'),'-0');
cfg.true_size               = round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
cfg.version                 = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_beta_%5.4f_R_gamma_%4.3f_dtri%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
                                cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.gamma,cfg.dTri,cfg.smooth,cfg.reinitR,cfg.add_source,cfg.dilate,cfg.niter_pcg);
cfg.out_dir                 = sprintf('./test_results/%s/%s',cfg.tag,cfg.version);

% broadcast dti_enhanced, there is not staggerred face-centered, only cell-centered
cfg.dti_enhanced = repmat(double(cfg.dti_enhanced), prod(cfg.true_size), 1);



for i = 1:global_steps
    cfg.vol(i).data = vol_tmp{i};
    rho_n = cfg.vol(i).data;
    if(~exist(sprintf('%s/rho_gt',cfg.out_dir), "dir"))
    mkdir(sprintf('%s/rho_gt',cfg.out_dir));
    end
    rho_gt_file = sprintf('%s/rho_gt/rhoGT_%s_t_%d.mat',cfg.out_dir,cfg.tag, cfg.first_time + i - 1);
    save(rho_gt_file, 'rho_n');
end

% --- Find global max signal for thresholding ---
data_max = zeros(size(cfg.vol(1).data));
for l = 1:length(cfg.vol)
    data_max = max(data_max, cfg.vol(l).data);
end

cfg.domain_size             = size(cfg.vol(1).data);

%% Run rOMT

if ~cfg.only_post_processing
    if cfg.reinitR == 1
        % --- Robust Parallel Pool Check ---
        % Check if a parallel pool is running, if not, start one.
        % This makes the script more robust for batch execution with nohup.
        pool = gcp('nocreate'); % Check for a pool without creating one
        if isempty(pool)
            fprintf('No parallel pool found. Starting a new one...\n');
            parpool; % If no pool, start one
        end        % --- End of Check ---
        [cfg, flag] = runROMT_par(cfg);
    else
        [cfg, flag] = runROMT(cfg);
    end
end

%% Try to load rOMT velocity results
fprintf('Loading rOMT results from: %s\n', cfg.out_dir);

global_steps = length(cfg.first_time:cfg.time_jump:cfg.last_time);
u_cell = cell(1, global_steps);

parfor tind = 1:global_steps
    ti = cfg.first_time+(tind-1)*cfg.time_jump;
    tf = cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump;
    ufile = sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,ti,tf,tind);

    if exist(ufile, 'file')
        fprintf('Loading %s\n', ufile);
        loaded_data = load(ufile, 'u');
        
        % --- MODIFICATION: Mask velocity based on density threshold ---
        rho = cfg.vol(tind).data;
        mask_low_density = rho < cfg.density_percent_thres;
        % WARNING: here abandon the first 4-7 time steps, because they do not preserve mass.
        if ismember(ti, cfg.exclude_frames)
            mask_low_density = true(size(mask_low_density)); % Exclude frames, set mask to true
        end
        % Reshape u for easier manipulation
        u_reshaped = reshape(loaded_data.u, [prod(cfg.domain_size), 3, cfg.nt]);
        
        % Apply mask to all velocity components and time steps
        u_reshaped(mask_low_density(:), :, :) = 0;
        
        u_cell{tind} = reshape(u_reshaped, size(loaded_data.u));
        % --- END MODIFICATION ---
    else
        error('Saved result file not found: %s\nCannot proceed to GLAD post-processing.', ufile);
    end
end
cfg.u = u_cell;

%% Run GLAD post-processing

[cfg, map, SL, stream, PATH] = runGLAD(cfg);

% [cfg, s, SL, PATH] = runGLAD2(cfg); % just faster analysis, without pecklet number.

% --- WARNING: could do Post-GLAD Pathline Filtering later ---
old_ind_brain = PATH.ind_brain; % Store old indices for reference



% --- Calculate axis limits from valid pathline start points for consistent visualization ---
fprintf('Calculating axis limits from valid pathline start points for consistent visualization...\n');
valid_start_points = PATH.startp(PATH.ind_brain, :);
if isempty(valid_start_points)
    warning('No valid pathlines found inside the brain mask. Using full volume for axis limits.');
    x_plot_lims = [1, cfg.true_size(2)];
    y_plot_lims = [1, cfg.true_size(1)];
    z_plot_lims = [1, cfg.true_size(3)];
else
    % The columns in PATH.startp are [y, x, z]
    % Plot's X-axis is dimension 2, Y-axis is dimension 1
    x_plot_lims = [min(valid_start_points(:, 2)), max(valid_start_points(:, 2))];
    y_plot_lims = [min(valid_start_points(:, 1)), max(valid_start_points(:, 1))];
    z_plot_lims = [min(valid_start_points(:, 3)), max(valid_start_points(:, 3))];
    
    % Add a small buffer to the limits to prevent data points from being on the edge
    x_buffer = (x_plot_lims(2) - x_plot_lims(1)) * 0.05;
    y_buffer = (y_plot_lims(2) - y_plot_lims(1)) * 0.05;
    z_buffer = (z_plot_lims(2) - z_plot_lims(1)) * 0.05;
    
    x_plot_lims = [x_plot_lims(1) - x_buffer, x_plot_lims(2) + x_buffer];
    y_plot_lims = [y_plot_lims(1) - y_buffer, y_plot_lims(2) + y_buffer];
    z_plot_lims = [z_plot_lims(1) - z_buffer, z_plot_lims(2) + z_buffer];
end
% --- End of axis calculation ---


%% Visualization of Pe map (full)
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);

x_slices = round(linspace(1, cfg.true_size(2), cfg.speedmap_slice));
y_slices = round(linspace(1, cfg.true_size(1), cfg.speedmap_slice));
z_slices = round(linspace(1, cfg.true_size(3), cfg.speedmap_slice));

figure,
hs=slice(y,x,z,map.Pe_full,x_slices,y_slices,z_slices); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
custom_alpha = [0, linspace(0.4, 0.6, 99)];
alpha('color'),alphamap(custom_alpha)

title(sprintf('%s: Pe Map',cfg.tag), 'Interpreter', 'none', 'FontSize', cfg.vis_font_size);
grid off, box off, axis image
xlabel('x-axis','FontSize',cfg.vis_font_size),ylabel('y-axis','FontSize',cfg.vis_font_size),zlabel('z-axis','FontSize',cfg.vis_font_size)


colormap(jet); % Apply the modified colormap

% Set a fixed color range for the Peclet number to avoid outliers washing
% out the map. A range is chosen because Pe=1 is the critical
% point where advection and diffusion are balanced. This range highlights
% the transition between diffusion-dominated (Pe<1) and 
% advection-dominated (Pe>1) regions.
clim([0, 800]);

axis tight; % Use 'axis tight' to preserve the aspect ratio
xlim(x_plot_lims);
ylim(y_plot_lims);
zlim(z_plot_lims);

set(gcf,'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');

view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')

% Add colorbar with custom tick at Pe = 1
cb = colorbar;
cb.Label.String = 'Peclet Number (Pe)';
cb.Label.FontSize = cfg.vis_font_size;

grid on,

saveas(gcf, sprintf('%s/%s/%s_LagPe_E%02d_%02d.png',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
savefig(gcf, sprintf('%s/%s/%s_LagPe_E%02d_%02d.fig',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%


%% Save and Visualize Velocity Field at Each Frame

if ~cfg.only_post_processing
    % --- First, find the global maximum speed across all frames for a consistent legend ---
    fprintf('Calculating global speed range for consistent visualization...\n');
    global_max_speed = 0;
    n = cfg.true_size;
    dim = 3;
    for tind = 1:length(cfg.u)
        u_interval_matrix = reshape(cfg.u{tind}, [prod(n) * dim, cfg.nt]);
        for k = 1:cfg.nt
            vel_frame = u_interval_matrix(:, k);
            u_xyz = reshape(vel_frame, [prod(n), dim]);
            speed_vec = sqrt(sum(u_xyz.^2, 2));
            current_max_speed = max(speed_vec);
            if current_max_speed > global_max_speed
                global_max_speed = current_max_speed;
            end
        end
    end
    % If no motion is detected, set a default non-zero max speed to avoid errors with clim
    if global_max_speed == 0
        global_max_speed = 1.0;
        fprintf('Global maximum speed was 0. Set to %.4f to avoid visualization errors.\n', global_max_speed);
    else
        fprintf('Global maximum speed is: %.4f\n', global_max_speed);
    end
    % --- Call the modular visualization function ---
    visualize_velocity_field(cfg, global_max_speed);
end

%% Visualization of speed map (at integral + avearage part)

% The speed map (map.s_full) represents the magnitude of the velocity at each voxel, averaged over time (integral part). 
% It's derived directly from the velocity field.
% Average from total step: 1 + (number of t1 steps) * nt * glacfg.nEulStep
s = map.s_full; % Speed map at integral part

figure,
hs=slice(y,x,z,s,x_slices,y_slices,z_slices); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap(linspace(0,1,100))
title(sprintf('%s: Speed Map',cfg.tag), 'Interpreter', 'none', 'FontSize', cfg.vis_font_size);
grid off, box off, axis image
xlabel('x-axis','FontSize',cfg.vis_font_size),ylabel('y-axis','FontSize',cfg.vis_font_size),zlabel('z-axis','FontSize',cfg.vis_font_size)
colormap(jet);

% Set a fixed color range for the speed map to [0, 1] to ensure
% consistency across different datasets and prevent high-velocity outliers
% from dominating the color scale.
clim([0, 0.8]);

axis tight; % Use 'axis tight' to preserve the aspect ratio
xlim(x_plot_lims);
ylim(y_plot_lims);
zlim(z_plot_lims);

set(gcf,'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');

view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
colorbar, grid on,

saveas(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
savefig(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.fig',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

fprintf("Speed map visualization saved.\n");

%% Visualization of flux vectors

% represent the displacement of particles over a certain time period. This displacement is calculated by integrating the velocity field along pathlines.
% quiver direction: the net direction of particle displacement
% quiver length: magnitude of the net displacement, also known as the pathline's displacement length
% quiver color: The color of each arrow also represents the magnitude of the net displacement (the same quantity as the length).

close all;
figure,

magnify = 1;%1;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

grid on;
axis equal;
axis tight;
hold on,

strid = cfg.strid; % Stride: plot one vector every 4 voxels
% quiver direction is the net direction of particle displacement
% quiver length is the magnitude of the net displacement from ti to tf
q = quiver3(PATH.startp(PATH.ind_brain(1:strid:end),2),PATH.startp(PATH.ind_brain(1:strid:end),1),PATH.startp(PATH.ind_brain(1:strid:end),3),PATH.disp(PATH.ind_brain(1:strid:end),2)*magnify, ...
    PATH.disp(PATH.ind_brain(1:strid:end),1)*magnify, ...
    PATH.disp(PATH.ind_brain(1:strid:end),3)*magnify, ... 
    'color','r','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','flux vectors');

title(sprintf('%s: velocity flux vectors',cfg.tag),'FontSize',cfg.vis_font_size, 'Interpreter', 'none'), set(gcf, 'Position', [376 49 1256 719])
xlabel('x-axis','FontSize',cfg.vis_font_size),ylabel('y-axis','FontSize',cfg.vis_font_size),zlabel('z-axis','FontSize',cfg.vis_font_size)
%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), q.WData(:)).^2, 2));

if numel(mags) == 0
    error('No valid pathlines were rendered, all pathlines are masked out ! please dilate the mask or turn the stride smaller.');
end

colormap_to_use = jet(256); % Use a standard-sized colormap

%// Get the current colormap
[~, ~, color_indices] = histcounts(mags, size(colormap_to_use, 1));

% cmap = jet(numel(q.UData)); % no used.
color_indices(color_indices == 0) = 1; % Handle edge cases where histcounts returns 0
arrow_colors = colormap_to_use(color_indices, :);

%  Create the correctly sized ColorData array for heads and tails
head_colors = repelem(arrow_colors, 3, 1); % 3 vertices per head
head_colors_rgba = uint8([head_colors, ones(size(head_colors,1), 1)]' * 255);

tail_colors = repelem(arrow_colors, 2, 1); % 2 vertices per tail
tail_colors_rgba = uint8([tail_colors, ones(size(tail_colors,1), 1)]' * 255);

set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', head_colors_rgba);
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', tail_colors_rgba);

%// Set the color of the quiver3 arrows
% legend('Location','best','Fontsize',12); % do not set because it is only mask + flux vector.
colorbar;
view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end

set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off');

clim([min(mags), max(mags)]);
set(gcf,'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');

saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%{
%% Visualization of advective & diffusive vectors
[InDadv,~,~] = intersect(PATH.ADVind,PATH.ind_brain);
[InDdiff,~,~] = intersect(PATH.DIFFind,PATH.ind_brain);

figure,
magnify = .5;%1;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

grid on, axis image
hold on,
q = quiver3(PATH.startp(InDadv(1:strid:end),2),PATH.startp(InDadv(1:strid:end),1),PATH.startp(InDadv(1:strid:end),3),PATH.disp(InDadv(1:strid:end),2)*magnify,PATH.disp(InDadv(1:strid:end),1)*magnify,PATH.disp(InDadv(1:strid:end),3)*magnify,'color','m','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','Advecitve Vectors');
title(sprintf('%s: Advective and Diffusive flux vectors',cfg.tag),'FontSize',20, 'Interpreter', 'none'), set(gcf, 'Position', [376 49 1256 719])
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)

hold on,
l = quiver3(PATH.startp(InDdiff(1:strid:end),2),PATH.startp(InDdiff(1:strid:end),1),PATH.startp(InDdiff(1:strid:end),3),PATH.disp(InDdiff(1:strid:end),2)*magnify,PATH.disp(InDdiff(1:strid:end),1)*magnify,PATH.disp(InDdiff(1:strid:end),3)*magnify,'color','g','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','Diffusive Vectors');

view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end
legend('Location','best','Fontsize',12);
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')

saveas(gcf, sprintf('%s/%s/%s_LagAdvDiffVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%}



fprintf("Flux vector visualization saved.\n");

%% Visualization of pathlines
figure,
SL2 = SL(PATH.ind_brain);
nSL = length(SL2);

strid = cfg.strid;

for ind = 1:strid:nSL
    SL_tmp = SL2{ind};
    colors = jet(size(SL_tmp,1));
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    %set(hlines,'EdgeColor',colors.round(PATH.displen(ind),:));
    set(hlines,'FaceVertexCData',[colors;colors(end,:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;

[x, y, z] = meshgrid(1:size(cfg.msk, 2), 1:size(cfg.msk, 1), 1:size(cfg.msk, 3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.2,.57,.2];
mskp.FaceAlpha= 0.05;
mskp.EdgeColor = [.2,.57,.2];
mskp.EdgeAlpha= 0;

if isfield(cfg, "msk_brain")
    [x, y, z] = meshgrid(1:size(cfg.msk_brain,2), 1:size(cfg.msk_brain,1), 1:size(cfg.msk_brain,3));
    mskfv = isosurface(x,y,z,cfg.msk_brain,0.5);
    mskp = patch(mskfv);
    mskp.FaceColor = [.57,.2,.2];
    mskp.FaceAlpha= 0.05;
    mskp.EdgeColor = [.57,.2,.2];
    mskp.EdgeAlpha= 0;
end

view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end
axis equal; % Use 'axis equal' to preserve the aspect ratio
axis tight;
set(gcf,'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');

colormap('jet'); colorbar, grid on,
title(sprintf('%s: Pathlines Starting with brain',cfg.tag), 'Interpreter', 'none', 'FontSize',cfg.vis_font_size);
saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump));


fprintf("Pathline visualization saved.\n");


% manually delete par pool, but if run again, do not close.
% delete(gcp('nocreate')) % delete existing parallel pool if any

