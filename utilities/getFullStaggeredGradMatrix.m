function G_full_staggered = getFullStaggeredGradMatrix(h1, h2, h3)
%GETFULLSTAGGEREDGRADMATRIX Creates a full gradient operator on a staggered grid.
%
%   G = getFullStaggeredGradMatrix(h1, h2, h3) creates a matrix G that,
%   when multiplied by a vector of scalar values at cell centers, produces
%   a vector of [gx; gy; gz] gradient components evaluated at *every* face
%   center (x-faces, then y-faces, then z-faces).
%
%   INPUTS:
%       h1, h2, h3: Vectors of cell sizes for x, y, z dimensions.
%
%   OUTPUT:
%       G: A sparse matrix that maps N cell-centered values to
%          3*Nfx + 3*Nfy + 3*Nfz face-centered gradient components.

    n1 = length(h1);
    n2 = length(h2);
    n3 = length(h3);

    % --- 1. Create Basic 1D Operators ---

    % 1a. 1D Nodal Differentiation (cell centers to faces)
    % This is the 'nodal' option from your ddx function.
    % n1 + 1 means rows, with first and last rows being ghost points.
    os = sparse(1,length(h));
    D_nodal_x = spdiags([-ones(n1,1), ones(n1,1)], [0, 1], n1 - 1, n1) ./ h1(1);
    D_nodal_x = [os; D_nodal_x; os];
    D_nodal_y = spdiags([-ones(n2,1), ones(n2,1)], [0, 1], n2 - 1, n2) ./ h2(1);
    D_nodal_y = [os; D_nodal_y; os];
    D_nodal_z = spdiags([-ones(n3,1), ones(n3,1)], [0, 1], n3 - 1, n3) ./ h3(1);
    D_nodal_z = [os; D_nodal_z; os];

    % 1b. 1D Cell-Centered Differentiation (cell centers to cell centers)
    % This is the 'ccn' option from your ddx function (internal part).
    D_center_x = spdiags([-ones(n1,1), ones(n1,1)], [-1, 1], n1, n1) ./ (2*h1(1));
    % because not in boundary, we could use simple forward/backward difference here
    D_center_x(1, 1:2)     = [-1 1] / h1(1);
    D_center_x(n1, n1-1:n1) = [-1 1] / h1(n1);
    
    D_center_y = spdiags([-ones(n2,1), ones(n2,1)], [-1, 1], n2, n2) ./ (2*h2(1));
    % because not in boundary, we could use simple forward/backward difference here
    D_center_y(1, 1:2)     = [-1 1] / h2(1);
    D_center_y(n2, n2-1:n2) = [-1 1] / h2(n2);

    D_center_z = spdiags([-ones(n3,1), ones(n3,1)], [-1, 1], n3, n3) ./ (2*h3(1));
    D_center_z(1, 1:2)     = [-1 1] / h3(1);
    D_center_z(n3, n3-1:n3) = [-1 1] / h3(n3);

    % 1c. 1D Averaging Operator (cell centers to faces)
    A_x = spdiags([ones(n1,1), ones(n1,1)], [0, 1], n1 + 1, n1) * 0.5;
    % ***Crucially, define the boundary conditions for A_x***
    % A common choice is extrapolation for the boundary faces:
    A_x(1, 1) = 1.0;         % Value at left-most x-face is same as u_center(1)
    A_x(end, n1) = 1.0;      % Value at right-most x-face is same as u_center(n1)

    A_y = spdiags([ones(n2,1), ones(n2,1)], [0, 1], n2 + 1, n2) * 0.5;
    A_y(1, 1) = 1.0;
    A_y(end, n2) = 1.0;
    
    A_z = spdiags([ones(n3,1), ones(n3,1)], [0, 1], n3 + 1, n3) * 0.5;
    A_z(1, 1) = 1.0;
    A_z(end, n3) = 1.0;

    % --- 2. Build 3D Operators using Kronecker Products ---
    % Helper for Kronecker product order
    kron3 = @(A,B,C) kron(A, kron(B,C));
    N_cells = n1 * n2 * n3;

    % --- 2a. Operators for Gradients at X-FACES ---
    Gxx = kron3(speye(n3), speye(n2), D_nodal_x);
    Gxy = kron3(speye(n3), D_center_y, A_x);
    Gxz = kron3(D_center_z, speye(n2), A_x);
    
    % Interleave the components for X-faces
    Nfx = size(Gxx, 1);
    G_at_X_faces = spalloc(3 * Nfx, N_cells, nnz(Gxx) + nnz(Gxy) + nnz(Gxz));
    G_at_X_faces(1:3:end, :) = Gxx; % gx components
    G_at_X_faces(2:3:end, :) = Gxy; % gy components
    G_at_X_faces(3:3:end, :) = Gxz; % gz components

    % --- 2b. Operators for Gradients at Y-FACES ---
    Gyx = kron3(speye(n3), A_y, D_center_x);
    Gyy = kron3(speye(n3), D_nodal_y, speye(n1));
    Gyz = kron3(D_center_z, A_y, speye(n1));

    % Interleave the components for Y-faces
    Nfy = size(Gyy, 1);
    G_at_Y_faces = spalloc(3 * Nfy, N_cells, nnz(Gyx) + nnz(Gyy) + nnz(Gyz));
    G_at_Y_faces(1:3:end, :) = Gyx; % gx components
    G_at_Y_faces(2:3:end, :) = Gyy; % gy components
    G_at_Y_faces(3:3:end, :) = Gyz; % gz components

    % --- 2c. Operators for Gradients at Z-FACES ---
    Gzx = kron3(A_z, speye(n2), D_center_x);
    Gzy = kron3(A_z, D_center_y, speye(n1));
    Gzz = kron3(D_nodal_z, speye(n2), speye(n1));

    % Interleave the components for Z-faces
    Nfz = size(Gzz, 1);
    G_at_Z_faces = spalloc(3 * Nfz, N_cells, nnz(Gzx) + nnz(Gzy) + nnz(Gzz));
    G_at_Z_faces(1:3:end, :) = Gzx; % gx components
    G_at_Z_faces(2:3:end, :) = Gzy; % gy components
    G_at_Z_faces(3:3:end, :) = Gzz; % gz components

    % --- 3. Stack all operators to form the final matrix ---
    G_full_staggered = [G_at_X_faces; G_at_Y_faces; G_at_Z_faces];

end