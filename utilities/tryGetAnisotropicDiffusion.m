function D_aniso = tryGetAnisotropicDiffusion(Grad, sigma, D_tensor, dti_enhanced)
%TRYGETANISOTROPICDIFFUSION_CENTERED Computes the anisotropic diffusion
%   operator using a cell-centered finite difference scheme.
%
%   INPUTS:
%       Grad:         A [3*N x N] cell-centered gradient matrix.
%       sigma:        A [N x 1] scalar field diffusion coefficient.
%       D_tensor:     An [N, 3, 3] stack of cell-centered tensors.
%       dti_enhanced: A scalar enhancement factor.
%
%   OUTPUT:
%       D_aniso: The final [N x N] sparse anisotropic diffusion operator.
% --- Handle the case of isotropic diffusion ---
if isempty(D_tensor)
    D_aniso = - sigma .* (Grad' * Grad);
    return;
end

N_voxels = size(D_tensor, 1);

% --- Construct the Block-Diagonal Tensor Operator (D_block) ---

% Scale the tensor stack by the scalar enhancement factor.
% D_tensor is already [N x 3 x 3], so this is a simple scalar multiplication.
D_enhanced_stack = dti_enhanced .* D_tensor; % Now [N_voxels x 3 x 3]

% This construction builds the sparse block-diagonal matrix efficiently.
V = D_enhanced_stack(:);

[row_template, col_template] = ndgrid(1:3, 1:3);

block_offsets = kron(3 * (0:N_voxels-1)', ones(9, 1));

I = repmat(row_template(:), N_voxels, 1) + block_offsets;
J = repmat(col_template(:), N_voxels, 1) + block_offsets;

% D_block is a [3*N x 3*N] sparse block-diagonal matrix
D_block = sparse(I, J, V, 3*N_voxels, 3*N_voxels);

D_aniso = - (Grad' * (D_block * Grad));

end