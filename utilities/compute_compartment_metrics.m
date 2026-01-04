function stats = compute_compartment_metrics(mask_opt, s_full, ds_full, Pe_full, SL, cfg, scount, do_print)
% in GLAD, Compute avg s_full, ds_full, Pe_full and v-flux inside a specific mask.
% mask_opt: struct with .nii path: 3D logical or numeric mask (same space as s_full/ds_full/Pe_full), threshold and name
% s_full, ds_full, Pe_full: 3D maps from GLAD (already averaged over visits)
% SL: cell array of pathlines in MATLAB grid coords [row(y), col(x), z]
% cfg: struct with fields true_size (n) and dx (mm)
% scount: optional visit count map (same size), improves "valid voxel" selection
% do_print: optional, default true

mask = nii2mat(mask_opt.path, cfg.x_range, cfg.y_range, cfg.z_range);
mask = mask > mask_opt.threshold;
if cfg.do_resize
    % Preserve binary mask on resize
    mask = resizeMatrix(single(mask), round(cfg.size_factor.*size(mask)), 'nearest') > 0.5;
end

if nargin < 8, do_print = true; end
n = size(s_full);

% Select valid voxels for averaging
if nargin >= 7 && ~isempty(scount)
    valid = mask(:) & (scount(:) > 0) & isfinite(s_full(:)) & (s_full(:) > 0);
else
    valid = mask(:) & isfinite(s_full(:)) & (s_full(:) > 0);
end

stats = struct();
stats.num_valid_vox = sum(valid);
if stats.num_valid_vox > 0
    stats.avg_s_full  = mean(s_full(valid));
    stats.avg_ds_full = mean(ds_full(valid));
else
    stats.avg_s_full = 0;
    stats.avg_ds_full = 0;
end

% Pe average (only where finite and > 0 within mask)
valid_pe = mask(:) & isfinite(Pe_full(:)) & (Pe_full(:) > 0);
stats.num_valid_pe_vox = sum(valid_pe);
if stats.num_valid_pe_vox > 0
    stats.avg_Pe_full = mean(Pe_full(valid_pe));
else
    stats.avg_Pe_full = 0;
end

% v-flux inside mask (unique visited voxels within mask)
visited_mask = false(n);
for i = 1:numel(SL)
    vxl = SL{i}; % [row(y), col(x), z] in grid units
    if isempty(vxl), continue; end
    vxl = vxl(all(isfinite(vxl),2),:); % drop NaNs/Infs
    vxl = round(vxl);

    % swap back correctly: first column is rows (y), second is cols (x)
    rows = max(1, min(vxl(:,1), n(1))); % y -> rows
    cols = max(1, min(vxl(:,2), n(2))); % x -> cols
    slcs = max(1, min(vxl(:,3), n(3))); % z -> slices

    lin  = sub2ind(n, rows, cols, slcs);
    visited_mask(lin) = true;
end
visited_in_mask = visited_mask & mask;

% Visualize three masks and block until closed
try
    fig = figure('Name', sprintf('Mask debug: %s', mask_opt.name), 'Color', [0.95 0.95 0.98]);
    hold on;
    [x, y, z] = meshgrid(1:n(2), 1:n(1), 1:n(3));

    if isfield(cfg, 'msk') && isequal(size(cfg.msk), n)
        mskfv = isosurface(x, y, z, cfg.msk, 0.5);
        if ~isempty(mskfv.vertices)
            p = patch(mskfv);
            p.FaceColor = [0.2, 0.57, 0.2]; p.FaceAlpha = 0.08;
            p.EdgeColor = [0.2, 0.57, 0.2]; p.EdgeAlpha = 0.0;
        end
    end

    % Compartment mask (from mask_opt)
    mskfv2 = isosurface(x, y, z, mask, 0.5);
    if ~isempty(mskfv2.vertices)
        p2 = patch(mskfv2);
        p2.FaceColor = [0.57, 0.2, 0.2]; p2.FaceAlpha = 0.25;
        p2.EdgeColor = [0.57, 0.2, 0.2]; p2.EdgeAlpha = 0.0;
    end

    % Plot visited_mask (not intersected) as requested
    visfv = isosurface(x, y, z, visited_in_mask, 0.5);
    if ~isempty(visfv.vertices)
        p3 = patch(visfv);
        p3.FaceColor = [0.2, 0.2, 0.57]; p3.FaceAlpha = 0.35;
        p3.EdgeColor = [0.2, 0.2, 0.57]; p3.EdgeAlpha = 0.0;
    end

    % Axes/view styling
    axis equal; axis tight; grid on;
    if isfield(cfg, 'view_azi_elevation'), view(cfg.view_azi_elevation); end
    if isfield(cfg, 'flip_z') && cfg.flip_z, set(gca, 'ZDir', 'reverse'); end
    title(sprintf('Masks: background (green), compartment (red), visited (blue)\n%s', mask_opt.name), 'Interpreter', 'none');

    waitfor(fig); % block until closed
catch ME
    warning(ME.identifier, '%s', ME.message);
end

stats.visited_voxel_count = nnz(visited_in_mask);
stats.voxel_volume_mm3    = cfg.dx^3;
stats.v_flux_mm3          = stats.visited_voxel_count * stats.voxel_volume_mm3;
stats.visited_fraction    = stats.visited_voxel_count / nnz(mask);

if do_print
    fprintf('mask=%s, valid=%d vox, avg speed=%.6f mm/min, avg diffuse_speed=%.6f mm/min, avg Pe=%.6f, v-flux=%.4f mm^3 (vox=%d)\n', ...
        mask_opt.name, stats.num_valid_vox, stats.avg_s_full * cfg.dx, stats.avg_ds_full * cfg.dx, stats.avg_Pe_full, stats.v_flux_mm3, stats.visited_voxel_count);
end
end