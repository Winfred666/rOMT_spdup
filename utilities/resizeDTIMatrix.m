function resized_exp_D = resizeDTIMatrix(D_tensor_5D, new_size, interp_mode)

% --- Vectorized Log-Euclidean Interpolation ---
disp('Using Log-Euclidean interpolation...');

% Get dimensions
[n1, n2, n3, ~, ~] = size(D_tensor_5D);

% Reshape for vectorized operations
D_flat = reshape(D_tensor_5D, n1*n2*n3, 3, 3);

% Symmetrize and handle bad tensors
D_flat = (D_flat + permute(D_flat, [1 3 2])) / 2;

% --- Part 1: Initial Cleanup ---

fro_norm = squeeze(sqrt(sum(D_flat.^2, [2,3])));
bad_indices = isnan(fro_norm) | isinf(fro_norm) | (fro_norm < 1e-12);
D_flat(bad_indices, :, :) = 0;
num_voxels = size(D_flat, 1);

% --- Part 2: Calculate the Matrix Logarithm (logm) ---
% pageeig/pagefun expect a different data layout, A for-loop is the clearest and most robust way to do this correctly.

log_D_flat = zeros(size(D_flat)); % Pre-allocate for speed

for i = 1:num_voxels
    % Get the 3x3 tensor for this voxel
    tensor = squeeze(D_flat(i, :, :));
    
    % Perform eigendecomposition on the single 3x3 matrix
    [V, Lambda] = eig(tensor, 'vector');
    
    % Prevent log of non-positive eigenvalues
    Lambda(Lambda <= 0) = 1e-12;
    
    % Reconstruct the log-tensor: V * log(Lambda_diag) * V'
    log_Lambda_diag = diag(log(Lambda));
    log_D_flat(i, :, :) = V * log_Lambda_diag * V';
end
log_D_flat(bad_indices, :, :) = 0; % Ensure bad tensors are zero in log space


% --- Part 3: Resize the Log-Tensor Field via Interpolation ---
% This part assumes you have the following variables defined:
%   - tensor_representation: The (n1 x n2 x n3 x 3 x 3) log-tensor field.
%   - new_size: A vector [new_n1, new_n2, new_n3] for the desired output size.
%   - interp_mode: (Optional) A string like '*linear' or '*cubic'.
tensor_representation = reshape(log_D_flat, n1, n2, n3, 3, 3);

% Set default interpolation mode if not provided
if ~exist('interp_mode', 'var')
    interp_mode = '*linear';
end

% Get original and new dimensions
[n1, n2, n3, ~, ~] = size(tensor_representation);
new_n1 = new_size(1);
new_n2 = new_size(2);
new_n3 = new_size(3);

disp(['Resizing log-tensor field from ' num2str(n1) 'x' num2str(n2) 'x' num2str(n3) ...
      ' to ' num2str(new_n1) 'x' num2str(new_n2) 'x' num2str(new_n3) '...']);

% Create normalized plaid grids for the original data points
% Handle singleton dimensions to avoid division by zero in normalization
x_orig = (n1 > 1) .* (0:n1-1)/(n1-1);
y_orig = (n2 > 1) .* (0:n2-1)/(n2-1);
z_orig = (n3 > 1) .* (0:n3-1)/(n3-1);
[x_mat, y_mat, z_mat] = ndgrid(x_orig, y_orig, z_orig);

% Create plaid grids for the new, desired interpolation points
x_new = (new_n1 > 1) .* (0:new_n1-1)/(new_n1-1);
y_new = (new_n2 > 1) .* (0:new_n2-1)/(new_n2-1);
z_new = (new_n3 > 1) .* (0:new_n3-1)/(new_n3-1);
[x_mat_interp, y_mat_interp, z_mat_interp] = ndgrid(x_new, y_new, z_new);

% Pre-allocate the output matrix for the resized log-tensors
resized_log_D = zeros([new_n1, new_n2, new_n3, 3, 3], 'like', tensor_representation);

% Interpolate each of the 6 unique tensor components (exploiting symmetry)
for i = 1:3
    for j = i:3
        % Extract the 3D volume for the current component (i,j)
        component_volume = tensor_representation(:, :, :, i, j);
        
        % Perform 3D interpolation on this component.
        % NOTE: interp3 expects grid vector arguments in (Y, X, Z) order.
        resized_component = interp3(y_mat, x_mat, z_mat, component_volume, ...
                                    y_mat_interp, x_mat_interp, z_mat_interp, ...
                                    interp_mode);
                                    
        % Handle potential NaNs from interpolation at the edges by setting them to 0
        resized_component(isnan(resized_component)) = 0;
        
        % Assign the resized component back, maintaining symmetry
        resized_log_D(:, :, :, i, j) = resized_component;
        resized_log_D(:, :, :, j, i) = resized_component;
    end
end

% The output 'resized_log_D' is now the resized tensor field in log-space.

% --- Part 4: Calculate the Matrix Exponential (expm) ---
% Again, a for-loop is the robust way to handle this for the concatenated faces.
num_voxel = size(resized_log_D, 1);
resized_exp_D = zeros(size(resized_log_D)); % Pre-allocate

for i = 1:num_voxel
    % Get the 3x3 log-tensor for this face
    log_tensor = squeeze(resized_log_D(i, :, :));
    % Perform eigendecomposition
    [V_voxel, L_voxel] = eig(log_tensor, 'vector');
    % Reconstruct the exp-tensor: V * exp(Lambda_diag) * V'
    exp_L_diag = diag(exp(L_voxel));
    resized_exp_D(i, :, :) = V_voxel * exp_L_diag * V_voxel';
end

% Optional: Save the resize tensors (No use now)