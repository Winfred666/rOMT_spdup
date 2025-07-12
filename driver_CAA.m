% This functions contains the rOMT algorithm and post-processing of the
% results on a sample dataset named 'C294' which is from a CAA dataset containsing about 48 rat brain cases w/w.o CAA disease.

addpath('./Inverse','./Sensitivities','./analyzeFlows',genpath('./utilities'))

%% set directories and parameters


cfg = set_config_CAA(); % Load configuration from function

% load ROI
if cfg.do_ROI_msk
    tmp = nii2mat(cfg.ROI_msk_path,cfg.x_range,cfg.y_range,cfg.z_range);
    cfg.msk = tmp>0;
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

<<<<<<< HEAD
% load vol
for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2
    tmp = nii2mat(sprintf('%s%02d%s',cfg.data_dir,cfg.first_time+(i-1)*cfg.time_jump,cfg.extension),cfg.x_range,cfg.y_range,cfg.z_range);
=======


% load vol

global_steps = (cfg.last_time-cfg.first_time)/cfg.time_jump+2;
vol_tmp = cell(global_steps, 1);

num_iter = global_steps;
parpool(num_iter); % Or parpool(12), etc.

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
>>>>>>> 9749637 (init)
    if cfg.do_resize
       tmp = resizeMatrix(tmp,round(cfg.size_factor.*size(tmp)),'linear');
    end
    if cfg.smooth>0
<<<<<<< HEAD
    tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    tmp(~cfg.msk) = 0;
    cfg.vol(i).data = tmp;
end

=======
        tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    tmp(~cfg.msk) = 0;
    vol_tmp{i} = tmp;
end



for i = 1:global_steps
    cfg.vol(i).data = vol_tmp{i};
end



>>>>>>> 9749637 (init)
cfg.domain_size             = size(cfg.vol(1).data);
cfg.sig_str                 = erase(num2str(cfg.sigma,'%.0e'),'-0');
cfg.true_size               = round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
cfg.version                 = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_beta_%5.4f_R_gamma_%4.3f_dtri%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
                                cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.gamma,cfg.dTri,cfg.smooth,cfg.reinitR,cfg.add_source,cfg.dilate,cfg.niter_pcg);
cfg.out_dir                 = sprintf('./test_results/%s/%s',cfg.tag,cfg.version);

%% Run rOMT
<<<<<<< HEAD
[cfg, flag] = runROMT(cfg);
%runROMT_par(cfg);

%% Run post-processing

%[cfg, map, SL, stream, PATH] = runGLAD(cfg);
[cfg, s, SL, PATH] = runGLAD2(cfg);


%% Visualization of speed map
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);
figure,
hs=slice(y,x,z,s,x,y,z); 
=======

if ~cfg.only_post_processing
    if cfg.reinitR == 1
        [cfg, flag] = runROMT_par(cfg);
    else
        [cfg, flag] = runROMT(cfg);
    end
end

% try to load previous results
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
        u_cell{tind} = loaded_data.u;
    else
        error('Saved result file not found: %s\nCannot proceed to GLAD post-processing.', ufile);
    end
end
cfg.u = u_cell;

%% Run GLAD post-processing

%[cfg, map, SL, stream, PATH] = runGLAD(cfg);

[cfg, s, SL, PATH] = runGLAD2(cfg);

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

%% Visualization of speed map (at integral part)
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);

x_slices = []; 
y_slices = []; 
z_slices = round(linspace(1, cfg.true_size(3), cfg.speedmap_slice));

figure,
hs=slice(y,x,z, s ,x_slices,y_slices,z_slices); 
>>>>>>> 9749637 (init)
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap(linspace(0,1,100))
title(sprintf('Test: tag = %s, Speed Map',cfg.tag),'Fontsize',20)
grid off, box off, axis image
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)
<<<<<<< HEAD
colormap(jet)
caxis([0,0.5])
view([-188.3500   13.7463])
=======
colormap(jet);

clim([0, 0.8*max(s(:))]);

axis equal; % Use 'axis equal' to preserve the aspect ratio
axis tight;
set(gcf,'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');

view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end
>>>>>>> 9749637 (init)
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
colorbar, grid on,

saveas(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
<<<<<<< HEAD
%% Visualization of flux vectors
figure,
strid = 10;
=======

fprintf("Speed map visualization saved.\n");

%% Visualization of flux vectors

% quiver direction: the net direction of particle displacement
% quiver length: magnitude of the net displacement, also known as the pathline's displacement length
% quiver color: The color of each arrow also represents the magnitude of the net displacement (the same quantity as the length).

close all;
figure,

>>>>>>> 9749637 (init)
magnify = .5;%1;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

<<<<<<< HEAD
grid on, axis image
hold on,
q = quiver3(PATH.startp(PATH.ind_brain(1:strid:end),2),PATH.startp(PATH.ind_brain(1:strid:end),1),PATH.startp(PATH.ind_brain(1:strid:end),3),PATH.disp(PATH.ind_brain(1:strid:end),2)*magnify,PATH.disp(PATH.ind_brain(1:strid:end),1)*magnify,PATH.disp(PATH.ind_brain(1:strid:end),3)*magnify,'color','r','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','flux vectors');
title(sprintf('%s: velocity flux vectors',cfg.tag),'FontSize',20, 'Interpreter', 'none'), set(gcf, 'Position', [376 49 1256 719])
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)
%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(jet);
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
legend('Location','best','Fontsize',12); colorbar;
view([-188.3500   13.7463]);
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
=======
grid on;
axis equal;
axis tight;
hold on,

strid = cfg.strid; % Stride: plot one vector every 4 voxels
q = quiver3(PATH.startp(PATH.ind_brain(1:strid:end),2),PATH.startp(PATH.ind_brain(1:strid:end),1),PATH.startp(PATH.ind_brain(1:strid:end),3),PATH.disp(PATH.ind_brain(1:strid:end),2)*magnify, ...
    PATH.disp(PATH.ind_brain(1:strid:end),1)*magnify, ...
    PATH.disp(PATH.ind_brain(1:strid:end),3)*magnify, ... 
    'color','r','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','flux vectors');

title(sprintf('%s: velocity flux vectors',cfg.tag),'FontSize',20, 'Interpreter', 'none'), set(gcf, 'Position', [376 49 1256 719])
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)
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
>>>>>>> 9749637 (init)

saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%{
%% Visualization of advective & diffusive vectors
[InDadv,~,~] = intersect(PATH.ADVind,PATH.ind_brain);
[InDdiff,~,~] = intersect(PATH.DIFFind,PATH.ind_brain);

figure,
<<<<<<< HEAD
strid = 10;
=======
>>>>>>> 9749637 (init)
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

<<<<<<< HEAD
view([-188.3500   13.7463]); legend('Location','best','Fontsize',12);
=======
view(cfg.view_azi_elevation);
if cfg.flip_z
    set(gca, 'ZDir', 'reverse');
end
legend('Location','best','Fontsize',12);
>>>>>>> 9749637 (init)
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')

saveas(gcf, sprintf('%s/%s/%s_LagAdvDiffVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%}
<<<<<<< HEAD
=======



fprintf("Flux vector visualization saved.\n");

>>>>>>> 9749637 (init)
%% Visualization of pathlines
figure,
SL2 = SL(PATH.ind_brain);
nSL = length(SL2);
<<<<<<< HEAD
%colors = jet(round(max(PATH.displen)));

for ind = 1:10:nSL
=======

strid = cfg.strid;

for ind = 1:strid:nSL
>>>>>>> 9749637 (init)
    SL_tmp = SL2{ind};
    colors = jet(size(SL_tmp,1));
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    %set(hlines,'EdgeColor',colors(round(PATH.displen(ind)),:));
    set(hlines,'FaceVertexCData',[colors;colors(end,:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;
<<<<<<< HEAD
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.37,.37,.37];
mskp.FaceAlpha= 0.1;
mskp.EdgeColor = [.37,.37,.37];
mskp.EdgeAlpha= 0;
view([-188.3500   13.7463])
axis image
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
colormap('jet'); colorbar, grid on,
title(sprintf('%s: Pathlines Starting with brain',cfg.tag), 'Interpreter', 'none', 'FontSize',18);
saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 









=======

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
title(sprintf('%s: Pathlines Starting with brain',cfg.tag), 'Interpreter', 'none', 'FontSize',18);
saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump));


fprintf("Pathline visualization saved.\n");

% manually delete par pool
delete(gcp('nocreate')) % delete existing parallel pool if any
>>>>>>> 9749637 (init)
