function cfg = set_config_CAA(config_path)
% set_config_CAA Load configuration for CAA experiments from a JSON file.
%   cfg = set_config_CAA(config_path)
%   - config_path: path to a JSON config under ./configs, or a base name
%                  (e.g., 'ours_ISO' -> './configs/ours_ISO.json').
%   The function applies sensible defaults, then overrides with user values.
%   dataset_name is forced to 'CAA'.

% Initialize defaults
cfg = struct();

% Always fixed
cfg.dataset_name            = 'CAA';

% Common defaults (can be overridden by JSON)
% cfg.only_post_processing    = 1;    % set to 1 if already run ROMT
% cfg.steady_velocity_file_path = '';% optional external steady velocity file

% Basic ROI/mask and volume defaults
cfg.do_ROI_msk              = 1;
cfg.ROI_msk_threshold       = 0.05;
cfg.x_range                 = 1:128;
cfg.y_range                 = 1:160;
cfg.z_range                 = 1:128;

% rOMT parameters (defaults similar to previous 'ours' case)
cfg.do_resize               = 1;
cfg.size_factor             = 0.6;

cfg.smooth                  = 1.8;      % diffusion smoothing time (multiple of dt)
cfg.dilate                  = 0;        % dilation for domain exploration
cfg.first_time              = 9;
cfg.time_jump               = 1;
cfg.last_time               = 26;

cfg.dti_enhanced            = 0.4;      % scale tensor field according to molecule size (1/3)
cfg.sigma                   = 2e-3;     % isotropic diffusion coefficient if no DTI

cfg.dx                      = 0.125;    % spatial resolution (mm)
cfg.dt                      = 1.0;      % time step (min)
cfg.nt                      = 4;        % sub-steps per time interval, 
% if MRI interval is 4min, make sure nt*dt = 4 for minute resolution 

cfg.gamma                   = 0.008;    % R_3 regularization (spatial smoothness)
cfg.beta                    = 0.0001;   % R_1 kinetic energy regularization
cfg.reinitR                 = 1;        % reinitialize rho each interval
cfg.reInitializeU           = 1;        % reinitialize u each interval
cfg.niter_pcg               = [550, 900];
cfg.maxUiter                = 20;
cfg.dTri                    = 1;        % boundary condition ('closed')
cfg.add_source              = 0;        % unbalanced rOMT disabled by default

% GLAD defaults
cfg.density_percent_thres   = 12;       % threshold to mask velocity by density
cfg.sp_thresh               = 12;       % start-point threshold (percent)
cfg.GLAD_spfs               = 18;       % sampling interval for start points
cfg.sl_tol                  = 2.5;      % minimum streamline length tol
cfg.GLAD_visualize_mask            = struct('name','brain_disabled','path','', 'threshold', 0.01);

cfg.integral_euler_step     = 10;       % Euler steps per dt for pathlines
cfg.GLAD_timestep_factor    = 1.0;      % extend pathline length factor

% Visualization defaults
cfg.speedmap_slice          = 12;
cfg.view_azi_elevation      = [5, 80];
cfg.strid                   = 1;
cfg.flip_z                  = 0;
cfg.vis_font_size           = 16;

% Normalize config_path input: allow bare name -> ./configs/<name>.json
if ~contains(config_path, filesep)
    candidate = fullfile('.', 'configs', [config_path, '.json']);
    if isfile(candidate)
        config_path = candidate;
    end
end
[~, cfg.tag, ~] = fileparts(config_path);

% If a directory was provided, error
if isfolder(config_path)
    error('set_config_CAA:InvalidPath', 'config_path must be a JSON file, not a directory: %s', config_path);
end

% If file extension missing but file not found, try appending .json
[config_dir, config_name, config_ext] = fileparts(config_path);
if isempty(config_ext)
    try_path = fullfile(config_dir, [config_name, '.json']);
    if isfile(try_path)
        config_path = try_path;
    end
end

% Load JSON and override defaults
if ~isfile(config_path)
    error('set_config_CAA:ConfigNotFound', 'Config file not found: %s', config_path);
end

raw = fileread(config_path);
user_cfg = jsondecode(raw);

% Override default fields with user-provided ones
uf = fieldnames(user_cfg);
for i = 1:numel(uf)
    f = uf{i};
    cfg.(f) = user_cfg.(f);
end

% Normalize range fields if provided in JSON
range_fields = {'x_range','y_range','z_range'};
for k = 1:numel(range_fields)
    fname = range_fields{k};
    if isfield(user_cfg, fname)
        cfg.(fname) = coerce_range(cfg.(fname));
    end
end

% Force dataset name
cfg.dataset_name = 'CAA';
cfg.anato                   = cfg.ROI_msk_path;

% Required fields validation
required_fields = {'data_template','ROI_msk_path'};
for i = 1:numel(required_fields)
    if ~isfield(cfg, required_fields{i}) || isempty(cfg.(required_fields{i}))
        error('set_config_CAA:MissingField', 'Required field missing in %s: %s', config_path, required_fields{i});
    end
end

% Optional: ensure exclude_frames exists
if ~isfield(cfg, 'exclude_frames')
    cfg.exclude_frames = [];
end

% Print brief confirmation
fprintf('Loaded CAA config: %s\n', config_path);

% Helper to coerce various JSON encodings into a MATLAB range vector
function r = coerce_range(v)
    if ischar(v) || (isstring(v) && isscalar(v))
        parts = strsplit(strtrim(string(v)), ':');
        nums = str2double(parts);
        if numel(nums) == 2
            r = nums(1):nums(2);
        elseif numel(nums) == 3
            r = nums(1):nums(2):nums(3);
        else
            error('Invalid range string: %s', v);
        end
    elseif isnumeric(v)
        v = double(v(:))';
        if numel(v) == 2
            r = v(1):v(2);
        elseif numel(v) == 3
            r = v(1):v(2):v(3);
        elseif numel(v) >= 1
            r = v; % explicit index list
        else
            error('Empty numeric range.');
        end
    elseif isstruct(v)
        if ~isfield(v,'start') || ~isfield(v,'stop')
            error('Range struct must have start and stop fields.');
        end
        step = 1;
        if isfield(v,'step') && ~isempty(v.step), step = v.step; end
        r = v.start:step:v.stop;
    else
        error('Unsupported range type for value: %s', class(v));
    end
    r = round(r);
    r = r(:)'; % row vector
end

end
