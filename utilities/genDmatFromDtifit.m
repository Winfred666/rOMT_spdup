% Add path to NIfTI toolbox if needed
addpath('./utilities');

base_path = '/data/xym/DTI_data/KX 078/dicom/out_dwi/dti_aligned';

% wasted original range in DWI space
% x_range = 1:128;
% y_range = 1:128;
% z_range = 1:31;

% Load eigenvectors (V1, V2, V3), should load using original method
tensor6 = load_untouch_nii(fullfile(base_path, 'dti_tensor.nii.gz')).img;

sz = size(tensor6);
if numel(sz) == 5
    fprintf("Reshape tensor6 of std SANS/ITK (X x Y x Z x 1 x 6) to D_tensor (XYZ x 3 x 3)\n");
    % Also from (X x Y x Z x 1 x 6) to (X x Y x Z x 3 x 3 for visualization)
    if sz(5) ~= 6
        error('tensor6 must be X x Y x Z x 1 x 6');
    end
    % WARNING: for MRtrix lower-triangle order
    Dxx = tensor6(:,:,:,:,1); % Dxx
    Dxy = tensor6(:,:,:,:,2); % Dxy
    Dyy = tensor6(:,:,:,:,3); % Dyy
    Dxz = tensor6(:,:,:,:,4); % Dxz
    Dyz = tensor6(:,:,:,:,5); % Dyz
    Dzz = tensor6(:,:,:,:,6); % Dzz
elseif numel(sz) == 4
    fprintf("Reshape tensor of dtifit X x Y x Z x 6 to D_tensor\n");
    if sz(4) ~= 6
        error('tensor6 must be X x Y x Z x 6');
    end
    % WARNING: for Dtifit upper-triangle order 
    Dxx = tensor6(:,:,:,1);
    Dxy = tensor6(:,:,:,2);
    Dxz = tensor6(:,:,:,3);
    Dyy = tensor6(:,:,:,4);
    Dyz = tensor6(:,:,:,5);
    Dzz = tensor6(:,:,:,6);
end


D_tensor_5D = zeros(sz(1), sz(2), sz(3), 3, 3, 'like', tensor6);
D_tensor_5D(:,:,:,1,1) = Dxx;
D_tensor_5D(:,:,:,2,2) = Dyy;
D_tensor_5D(:,:,:,3,3) = Dzz;
D_tensor_5D(:,:,:,1,2) = Dxy;
D_tensor_5D(:,:,:,2,1) = Dxy;
D_tensor_5D(:,:,:,1,3) = Dxz;
D_tensor_5D(:,:,:,3,1) = Dxz;
D_tensor_5D(:,:,:,2,3) = Dyz;
D_tensor_5D(:,:,:,3,2) = Dyz;

D_tensor = D_tensor_5D;
% PIC , cell centered , not staggered.
save(fullfile(base_path, 'dti_tensor_3_3.mat'), 'D_tensor');

% Fill 3x3 tensor for each voxel
D_tensor = reshape(D_tensor_5D, [], 3, 3); % (XYZ) x 3 x 3
save(fullfile(base_path, 'dti_tensor.mat'), 'D_tensor');
