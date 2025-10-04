function stats = compute_compartment_metrics(mask_opt, s_full, ds_full, Pe_full, SL, cfg, scount, do_print)
% in GLAD, Compute avg s_full, ds_full, Pe_full and v-flux inside a specific mask.
% mask_opt: struct with .nii path: 3D logical or numeric mask (same space as s_full/ds_full/Pe_full), threshold and name
% s_full, ds_full, Pe_full: 3D maps from GLAD (already averaged over visits)
% SL: cell array of pathlines in MATLAB grid coords [x(col), y(row), z], as in runGLAD
% cfg: struct with fields true_size (n) and dx (mm)
% scount: optional visit count map (same size), improves "valid voxel" selection
% do_print: optional, default true

mask = nii2mat(mask_opt.path, cfg.x_range, cfg.y_range, cfg.z_range);
mask = mask > mask_opt.threshold;
if cfg.do_resize
    mask = resizeMatrix(double(mask),round(cfg.size_factor.*size(mask)),'linear');
    mask(mask~=1) = 0; % make sure it is binary
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
    vxl = round(SL{i}); % [x, y, z] in grid units
    if isempty(vxl), continue; end
    vxl = unique(vxl, 'rows');
    cols = max(1, min(vxl(:,1), n(2))); % x -> columns
    rows = max(1, min(vxl(:,2), n(1))); % y -> rows
    slcs = max(1, min(vxl(:,3), n(3))); % z -> slices
    lin  = sub2ind(n, rows, cols, slcs);
    visited_mask(lin) = true;
end
visited_in_mask = visited_mask & mask;
stats.visited_voxel_count = nnz(visited_in_mask);
stats.voxel_volume_mm3    = cfg.dx^3;
stats.v_flux_mm3          = stats.visited_voxel_count * stats.voxel_volume_mm3;
stats.visited_fraction    = stats.visited_voxel_count / nnz(mask);

if do_print
    fprintf('mask=%s, valid=%d vox, avg speed=%.6f mm/min, avg diffuse_speed=%.6f mm/min, avg Pe=%.6f, v-flux=%.4f mm^3 (vox=%d)\n', ...
        mask_opt.name, stats.num_valid_vox, stats.avg_s_full * cfg.dx, stats.avg_ds_full * cfg.dx, stats.avg_Pe_full, stats.v_flux_mm3, stats.visited_voxel_count);
end
end