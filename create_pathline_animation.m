%% Pathline Animation Script
% This script iterates through a range of sp_thresh values, generates a
% pathline visualization for each, and compiles the resulting images into
% a GIF animation.
%
% Before running, ensure that the required variables (cfg, PATH, SL) are
% loaded in the MATLAB workspace from a previous analysis run.

% --- Configuration ---
sp_thresh_range = 4:50; % Range of thresholds to iterate through
output_gif_filename = sprintf('%s/%s/Pathline_Animation_E%02d_%02d.gif', cfg.out_dir, cfg.outdir_v, cfg.first_time, cfg.last_time+cfg.time_jump);
frame_delay = 0.1; % Delay between frames in the GIF (in seconds)

% --- Pre-calculation (for efficiency) ---
fprintf('Pre-calculating intensity data for animation...\n');
% Calculate data_max, the maximum intensity projection over time
data_max = zeros(size(cfg.vol(1).data));
for l = 1:length(cfg.vol)
    data_max = max(data_max, cfg.vol(l).data);
end

% Get the starting points of all pathlines and their intensities
brain_start_points = PATH.startp(:, :); % [y, x, z]
brain_start_points_rounded = round(brain_start_points);
brain_start_points_rounded = max(brain_start_points_rounded, 1);
dims = size(data_max);
brain_start_points_rounded(:,1) = min(brain_start_points_rounded(:,1), dims(1));
brain_start_points_rounded(:,2) = min(brain_start_points_rounded(:,2), dims(2));
brain_start_points_rounded(:,3) = min(brain_start_points_rounded(:,3), dims(3));
start_indices_linear = sub2ind(dims, brain_start_points_rounded(:,1), brain_start_points_rounded(:,2), brain_start_points_rounded(:,3));
start_point_intensities = data_max(start_indices_linear);

% --- Loop and Generate Frames ---
fprintf('Generating frames for animation...\n');
fig = figure('Visible', 'off'); % Create a figure but keep it hidden for speed

for i = 1:length(sp_thresh_range)
    sp_thresh = sp_thresh_range(i);
    fprintf('Processing sp_thresh = %d...\n', sp_thresh);
    
    % Filter pathlines based on the current intensity threshold
    passed_thresh_mask = find((start_point_intensities > sp_thresh) & (cfg.msk_brain(start_indices_linear) > 0));
    
    % --- Visualization ---
    clf(fig); % Clear the current figure for the new plot
    set(fig, 'unit','normalized','position',[0.1, 0.1, 0.6, 0.7],'Color',[0.85,0.85,0.93], 'InvertHardcopy', 'off');
    ax = axes(fig);
    hold(ax, 'on');
    
    SL2 = SL(passed_thresh_mask);
    nSL = length(SL2);
    strid = cfg.strid;

    for ind = 1:strid:nSL
        SL_tmp = SL2{ind};
        colors = jet(size(SL_tmp,1));
        hlines = patch(ax, [SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
        set(hlines,'FaceVertexCData',[colors;colors(end,:)],'EdgeColor','flat','FaceColor','none');
    end

    [x, y, z] = meshgrid(1:size(cfg.msk, 2), 1:size(cfg.msk, 1), 1:size(cfg.msk, 3));
    mskfv = isosurface(x,y,z,cfg.msk,0.5);
    mskp = patch(ax, mskfv);
    mskp.FaceColor = [.2,.57,.2]; mskp.FaceAlpha= 0.05; mskp.EdgeColor = 'none';

    if isfield(cfg, "msk_brain")
        [x, y, z] = meshgrid(1:size(cfg.msk_brain,2), 1:size(cfg.msk_brain,1), 1:size(cfg.msk_brain,3));
        mskfv = isosurface(x,y,z,cfg.msk_brain,0.5);
        mskp = patch(ax, mskfv);
        mskp.FaceColor = [.57,.2,.2]; mskp.FaceAlpha= 0.05; mskp.EdgeColor = 'none';
    end

    view(ax, cfg.view_azi_elevation);
    if cfg.flip_z
        set(ax, 'ZDir', 'reverse');
    end
    axis(ax, 'equal', 'tight');
    grid(ax, 'on');
    colormap(ax, 'jet');
    colorbar(ax);
    title(ax, sprintf('%s: Pathlines (sp\\_thresh = %d)',cfg.tag, sp_thresh), 'Interpreter', 'none', 'FontSize',cfg.vis_font_size);
    hold(ax, 'off');
    
    % Capture the frame for the GIF
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF File
    if i == 1
        imwrite(imind, cm, output_gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay);
    else
        imwrite(imind, cm, output_gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end

close(fig); % Close the hidden figure
fprintf('Animation saved to: %s\n', output_gif_filename);
