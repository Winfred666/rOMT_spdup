% Add path to NIfTI toolbox if needed
addpath('./utilities/NIfTI_analyze');

base_path = '/data/xym/DTI_data/ISO 052/dicom/out_dwi';
x_range = 1:128;
y_range = 1:128;
z_range = 1:31;
% Load eigenvalues
L1 = nii2mat(fullfile(base_path, 'dti_L1.nii.gz'), x_range, y_range, z_range);
L2 = nii2mat(fullfile(base_path, 'dti_L2.nii.gz'));
L3 = nii2mat(fullfile(base_path, 'dti_L3.nii.gz'));

% Load eigenvectors (V1, V2, V3)
V1 = nii2mat(fullfile(base_path, 'dti_V1.nii.gz'));
V2 = nii2mat(fullfile(base_path, 'dti_V2.nii.gz'));
V3 = nii2mat(fullfile(base_path, 'dti_V3.nii.gz'));

% Get image size
[dimX, dimY, dimZ, ~] = size(V1);

% Preallocate tensor array
tensor = zeros(dimX, dimY, dimZ, 3, 3);

% Loop over voxels
for x = 1:dimX
    for y = 1:dimY
        for z = 1:dimZ
            % Construct eigenvector matrix V (3x3)
            V = [squeeze(V1(x,y,z,:)), squeeze(V2(x,y,z,:)), squeeze(V3(x,y,z,:))];
            
            % Construct diagonal eigenvalue matrix L
            L = diag([L1(x,y,z), L2(x,y,z), L3(x,y,z)]);
            
            % Reconstruct tensor: D = V * L * V'
            tensor(x,y,z,:,:) = V * L * V';
        end
    end
end

D_tensor = tensor;

save(fullfile(base_path,'DTI_tensor_3_3.mat'), 'D_tensor', '-v7.3');

[dimX, dimY, dimZ, ~, ~] = size(tensor);
N = dimX * dimY * dimZ;
D_tensor = reshape(tensor, [N, 3, 3]);  % reshape to N x 3 x 3

save(fullfile(base_path,'DTI_tensor.mat'), 'D_tensor', '-v7.3');
