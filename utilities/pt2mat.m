function im_out = pt2mat(data_path,x1_range,x2_range,x3_range,size_factor,up_thresh,current_frame)
% this function extracts concentration data directly from a PyTorch .pt file.
% No clamping, smoothing, resize, or metadata parsing is performed.

narginchk(4, 7);
if nargin < 7
    current_frame = [];
end
unused_legacy_args = {size_factor, up_thresh}; %#ok<NASGU>

if ~isfile(data_path)
	error('pt2mat:FileNotFound', 'PT file not found: %s', data_path);
end

% Load tensor using the external Python executable specified in ROMT_PYTHON_EXE.
py_exe = '/usr/bin/python3.10';

if isempty(py_exe)
    error('pt2mat:NoPythonExe', 'Environment variable ROMT_PYTHON_EXE must be set to the full python executable path used by rOMT.');
end

try
    im_out = load_selected_volume_via_external_python(data_path, x1_range, x2_range, x3_range, current_frame, py_exe);
catch ME
    error('pt2mat:LoadFailed', ['Failed to load .pt file using external Python (%s): %s'], py_exe, ME.message);
end

end

function im_out = load_selected_volume_via_external_python(data_path, x1_range, x2_range, x3_range, current_frame, py_exe)
this_dir = fileparts(mfilename('fullpath'));
script_path = fullfile(this_dir, 'pt_extract_tensor.py');
if ~isfile(script_path)
	error('pt2mat:MissingHelperScript', 'Required helper script not found: %s', script_path);
end

tmp_root = tempname;
out_bin = [tmp_root, '.bin'];

x1_str = char(join(string(x1_range(:)'), ','));
x2_str = char(join(string(x2_range(:)'), ','));
x3_str = char(join(string(x3_range(:)'), ','));

cmd = sprintf('"%s" "%s" --input "%s" --out-bin "%s" --x1 "%s" --x2 "%s" --x3 "%s"', ...
	char(py_exe), char(script_path), char(data_path), char(out_bin), x1_str, x2_str, x3_str);
if ~isempty(current_frame)
	cmd = sprintf('%s --frame %d', cmd, current_frame);
end
[status, output] = system(cmd);
if status ~= 0
	error('pt2mat:ExternalPythonFailed', 'Command failed (status=%d): %s', status, output);
end

if ~isfile(out_bin)
	error('pt2mat:ExtractorOutputMissing', 'Extractor did not produce expected output: %s', out_bin);
end

cleanup_obj = onCleanup(@() cleanup_temp_outputs(out_bin)); %#ok<NASGU>

fid = fopen(out_bin, 'rb');
if fid < 0
	error('pt2mat:OpenBinFailed', 'Cannot open extractor output bin: %s', out_bin);
end
file_cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

numel_expected = numel(x1_range) * numel(x2_range) * numel(x3_range);
data = fread(fid, numel_expected, 'double=>double');
if numel(data) ~= numel_expected
	error('pt2mat:BinSizeMismatch', 'Expected %d doubles in %s, got %d.', numel_expected, out_bin, numel(data));
end

% Binary is C-order contiguous for shape [len(x1), len(x2), len(x3)].
im_out = reshape(data, [numel(x3_range), numel(x2_range), numel(x1_range)]);
im_out = permute(im_out, [3, 2, 1]);
im_out = double(im_out);
end

function cleanup_temp_outputs(out_bin)
if isfile(out_bin), delete(out_bin); end
end

