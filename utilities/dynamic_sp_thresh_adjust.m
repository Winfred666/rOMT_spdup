% Copy this code after running GLAD, to quickly getting pathline graph at a higher sp_thresh

sp_thresh = 24;

% --- Post-GLAD Pathline Filtering ---
fprintf('Performing post-GLAD pathline filtering based on sp_thresh...\n');
% Calculate data_max, which is the maximum intensity projection over time
data_max = zeros(size(cfg.vol(1).data));
for l = 1:length(cfg.vol)
    data_max = max(data_max, cfg.vol(l).data);
end

% Get the starting points of pathlines within the brain
brain_start_points = PATH.startp(:, :); % [y, x, z]

% Convert start points to linear indices
brain_start_points_rounded = round(brain_start_points);
brain_start_points_rounded = max(brain_start_points_rounded, 1); % Ensure indices are at least 1
dims = size(data_max);
brain_start_points_rounded(:,1) = min(brain_start_points_rounded(:,1), dims(1));
brain_start_points_rounded(:,2) = min(brain_start_points_rounded(:,2), dims(2));
brain_start_points_rounded(:,3) = min(brain_start_points_rounded(:,3), dims(3));

start_indices_linear = sub2ind(dims, brain_start_points_rounded(:,1), brain_start_points_rounded(:,2), brain_start_points_rounded(:,3));

% Get the intensity at each starting point
start_point_intensities = data_max(start_indices_linear);

% Filter pathlines based on the intensity threshold

passed_thresh_mask = find((start_point_intensities > sp_thresh) & (cfg.msk_brain(start_indices_linear) > 0));

PATH.ind_brain = passed_thresh_mask; % Update PATH.ind_brain with filtered indices




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
title(sprintf('%s: Pathlines Starting with brain thresh %02d',cfg.tag, sp_thresh), 'Interpreter', 'none', 'FontSize',cfg.vis_font_size);
saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump, sp_thresh));


fprintf("Pathline visualization saved.\n");

% manually delete par pool, but if run again, do not close.
% delete(gcp('nocreate')) % delete existing parallel pool if any
