function G = getCellCenteredGradMatrix(bcFlag, h1, h2, h3)
%GETCENTEREDGRADMATRIX Creates a cell-centered gradient operator.
%   G = getCenteredGradMatrix(h1, h2, h3) creates a matrix G that, when
%   multiplied by a vector of scalar values at cell centers, produces a
%   vector of [gx; gy; gz] gradient components, also at the cell centers.
%   INPUTS:
%       h1, h2, h3: Vectors of cell sizes for x, y, z dimensions.
%
%   OUTPUT:
%       G: A sparse [3*N x N] matrix, where N is the total number of cells.
% bcFlag could only be neumann boundary.
n1 = length(h1);
n2 = length(h2);
% --- Create 1D Central Difference Operators ---
% For a vector u, D*u computes du/dx at each point.
% We use forward/backward differences at the boundaries (Neumann BC).
% Dx operator for the x-direction
e = ones(n1, 1);
Dx_1D = spdiags([-e e], [-1 1], n1, n1) / (2*h1(1)); % Assuming constant h1
Dx_1D(1, 1:2) = [-1 1] / h1(1);   % Forward difference at first point
Dx_1D(n1, n1-1:n1) = [-1 1] / h1(n1); % Backward difference at last point

% Dy operator for the y-direction
e = ones(n2, 1);
Dy_1D = spdiags([-e e], [-1 1], n2, n2) / (2*h2(1)); % Assuming constant h2
Dy_1D(1, 1:2) = [-1 1] / h2(1);
Dy_1D(n2, n2-1:n2) = [-1 1] / h2(n2);
if nargin < 3
    error('You must supply the boundary condition flag: ''ccn'' or ''ccd'' for Cell Centered Nuemann and Cell Centered Dirchlet, respectively.')
elseif nargin == 3
    Gx = kron(speye(n2), Dx_1D);
    Gy = kron(Dy_1D, speye(n1));
    G = [Gx; Gy];
elseif nargin == 4
    n3 = length(h3);
    % Dz operator for the z-direction
    e = ones(n3, 1);
    Dz_1D = spdiags([-e e], [-1 1], n3, n3) / (2*h3(1)); % Assuming constant h3
    Dz_1D(1, 1:2) = [-1 1] / h3(1);
    Dz_1D(n3, n3-1:n3) = [-1 1] / h3(n3);

    % --- Extend to 3D using Kronecker products ---
    Gx = kron(speye(n3), kron(speye(n2), Dx_1D));
    Gy = kron(speye(n3), kron(Dy_1D, speye(n1)));
    Gz = kron(Dz_1D, kron(speye(n2), speye(n1)));

    % --- Interleave operators to create the full gradient matrix ---
    % This orders the output as [gx(1);gy(1);gz(1); gx(2);gy(2);gz(2); ...]
    % which is required for multiplication with a block-diagonal D_block.
    N_cells = n1 * n2 * n3;
    G = spalloc(3 * N_cells, N_cells, nnz(Gx) + nnz(Gy) + nnz(Gz));
    G(1:3:end, :) = Gx; % Place all gx rows in 1, 4, 7, ...
    G(2:3:end, :) = Gy; % Place all gy rows in 2, 5, 8, ...
    G(3:3:end, :) = Gz; % Place all gz rows in 3, 6, 9, ...

end

% %% Get Cell Centered Gradient Matrix
% %
% % Returns the gradient matix in either two or three dimensions.
% %
% %       G = getCellCenteredGradMatrix(bcFlag,h1,h2,h3);
% %
% % Where G is the gradient and h# = [d1 d2 ... dn]'; where d is the size of 
% % each cell.
% %
% % Here, 'bcFlag' is a boundary condition flag. Set this to either: 
% %
% %           'ccn'           for Cell Centered Nuemann
% %           'ccd'           for Cell Centered Dirchlet
% %
% % If bcFlag is a cell, then it is assumed that each element in the cell is
% % a string that sets the boundary conditions of that dimension of the
% % matrix.
% %
% % For example:
% %       
% %       bcFlag = {'ccn','ccd'}
% %       G = getCellCenteredGradMatrix(bcFlag,h1,h2);
% %
% % To get greater control over the boundary conditions use:
% %
% %   [G B] = getCellCenteredGradMatrix(bc,h1,h2,h3);
% %
% % Where B is the boundary conditions matrix, that puts values from a vector
% % into the same space as the gradient so you can add them to get the true
% % gradient. Note that because BC are know, you can actually move B over to
% % the RHS in practice.
% %
% % Rowan Cockett & Eldad Haber
% % 24-May-2012
% % University of British Columbia
% % rcockett@eos.ubc.ca
% %
% % See Also ddx testOperators testGradBC

% function [G,B,Gx,Gy,Gz] = getCellCenteredGradMatrix(bcFlag,h1,h2,h3)
% if ischar(bcFlag)
%     BC = {bcFlag,bcFlag,bcFlag};
% elseif iscell(bcFlag) && length(bcFlag) == (nargin-1)
%     BC = bcFlag;
% else
%     error('Different boundary conditions are supported, use a cell array. E.g. {''ccn'',''ccd''}')
% end
% if nargin < 3
%     error('You must supply the boundary condition flag: ''ccn'' or ''ccd'' for Cell Centered Nuemann and Cell Centered Dirchlet, respectively.')
% elseif nargin == 3
%     % Create the Gradient in 2D
%     n1  = length(h1); n2 = length(h2);
%     Gx = kron(speye(n2),ddx(h1(:),BC{1}));
%     Gy = kron(ddx(h2(:),BC{2}),speye(n1));
%     G  = [Gx;Gy];
%     if nargout > 1
%         Bx = kron(speye(n2),ddx(h1(:),[BC{1},'BC']));
%         By = kron(ddx(h2(:),[BC{2},'BC']),speye(n1));
%         B  = blkdiag(Bx, By);
%     end
% elseif nargin == 4
%     % Create the stagger Gradient in 3D !!!
%     % WARNING: to better aligned with DTI tensor field, need PIC Gradient.
%     n1  = length(h1); n2 = length(h2); n3 = length(h3);
%     kron3 = @(A,B,C) kron(A,kron(B,C));
%     Gx = kron3(speye(n3),speye(n2),ddx(h1(:),BC{1}));
%     Gy = kron3(speye(n3),ddx(h2(:),BC{2}),speye(n1));
%     Gz = kron3(ddx(h3(:),BC{3}),speye(n2),speye(n1));
%     G  = [Gx;Gy;Gz];
    
%     if nargout > 1
%         Bx = kron3(speye(n3),speye(n2),ddx(h1(:),[BC{1},'BC']));
%         By = kron3(speye(n3),ddx(h2(:),[BC{2},'BC']),speye(n1));
%         Bz = kron3(ddx(h3(:),[BC{3},'BC']),speye(n2),speye(n1));
%         B  = blkdiag(Bx, By, Bz);
%     end
% else
%     error('Only 2D and 3D supported, use ddx for 1D.');
% end



