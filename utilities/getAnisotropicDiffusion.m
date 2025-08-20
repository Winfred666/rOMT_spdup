function D_aniso = getAnisotropicDiffusion(Grad, D_tensor)
    % Grad: (3N x N), column stacked [Gx; Gy; Gz] for every voxel.
    % D_tensor: (N x 3 x 3) rows [Dxx Dxy Dxz; Dyx Dyy Dyz; Dzx Dzy Dzz]
    N = size(D_tensor,1);
    M = size(Grad,2);          % should equal N
    D_aniso = sparse(M,M);
    for i = 1:N
        rows = (i-1)*3 + (1:3); % 3 rows for each voxel
        Gi = Grad(rows,:);     % 3 x N
        % Di is the diffusion tensor for voxel i, 3x3 matrix
        Di = D_tensor(i,:,:); % 1 x 3 x 3
        Di = squeeze(Di);     % 3 x 3 matrix
        % Compute the contribution to the anisotropic diffusion matrix
        % This results in a contribution of size N x N
        D_aniso = D_aniso + Gi' * Di * Gi;
    end
end