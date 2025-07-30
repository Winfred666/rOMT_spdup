function [u,phi,g] = GNblock_u(Rho_i,u,nt,dt,par,tag_str)
%%
%% Note: Changes made by Xinan Chen in Apr 21: (modified 3D version in this project too)
% (1) in GNblock_u.m: move "U = reshape(u,2*prod(par.n),[]);" and "dmk2 = zeros(2*prod(par.n),nt);" into the loop
% (2) in get_dRudu.m: remove negative sign computing G
% (3) in get_drNduT.m: remove unnecessary variable sensi to avoid computing pcg twice
% (4) add new function get_drNduTdrNdu.m to combine
% get_drNduT(Rho_i,u,nt,dt,par,get_drNdu(Rho_i,u,nt,dt,par,x)) togeter in
% computing H
% (5) add new version in dTrilinear.m with equivalent calculation in math,
% even though no obvious running improvement has shown
% (6) add Rho as input in get_drNduTdrNdu.m and get_drNduT2.m, to avoid repeated computing of
% Rho = advecDiff.m
% (7) add new function get_drNduTdrNdu2.m and get_drNduT3.m to store S,Tx,Ty,Tz
% (8) combine 2d and 3d together. Modified functions include: GNblock_u.m,
% paramInitFunc.m, advecDiff.m, dTrilinear2d.m (new),  dTrilinear3d.m
% (new), get_drNduT3.m, get_drNduTdrNdu2.m
% (9) add new input dTri in paramInitFunc.m to pass type of boundary condition to
% dTrilinears2d.m and dTrilinears3d.m
% (10) in GNblock_u.m, precompute part of g and H
% (11) in paramInitFunc.m, precompute par.B = I-dt*Mdis and save to par,
% resulting in modifications in get_drNduT3.m, get_drNduTdrNdu2.m and
% advecDiff.m
% (12) create get_drNduTdrNdu3.m and get_drNduT4.m to further reduce
% running time
% (13) modified sdiag.m to slightly reduce running time

if nargin < 6
    tag_str = '';
end

% --- MODIFICATION: Initialize history arrays ---
loss_history.mk = [];
loss_history.phiN = [];
loss_history.Ru = [];
loss_history.total = [];
% --- END MODIFICATION ---

[phi, ~, ~, ~] = get_phi(Rho_i,u,par); % Use ~ to ignore unused outputs
A            = kron(ones(1,par.dim),speye(prod(par.n)));
Abig         = kron(speye(nt),A);
flag         = 0;

%% Loop
for i = 1:par.maxUiter
    Rho      = advecDiff(Rho_i,u,nt,dt,par);
    U        = reshape(u,par.dim*prod(par.n),[]);
    dmk2     = zeros(par.dim*prod(par.n),nt);
    
    RHO0 = [Rho_i,Rho(:,1:end-1)];
    if par.dim==2
        for k = 1:nt
            U1 = reshape(U(1:prod(par.n),k),par.n');
            U2 = reshape(U(prod(par.n)+1:end,k),par.n');
            [M.S{k},M.Tx{k},M.Ty{k}]  = dTrilinears2d(RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2,...
                             par.h1(1),par.h2(1),par.bc);
        end
    elseif par.dim==3
        for k = 1:nt
            U1 = reshape(U(1:prod(par.n),k),par.n');
            U2 = reshape(U(prod(par.n)+1:2*prod(par.n),k),par.n');
            U3 = reshape(U(2*prod(par.n)+1:end,k),par.n');
            [M.S{k},M.Tx{k},M.Ty{k},M.Tz{k}]  = dTrilinears3d(RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                             par.h1(1),par.h2(1),par.h3(1),par.bc);
        end
    else
        warning('In GNblock_u.m: dimension of data should be either 2 or 3')
    end   
    
    for      j = 1:nt
             dmk2(:,1:j) = dmk2(:,1:j) + reshape(par.hd*dt*get_drNduT4(M,j,dt,par,A*(U(:,j).*U(:,j))),par.dim*prod(par.n),j);
    end
    
    g        = (par.beta*2*par.hd*dt*Rho(:)'*Abig*sdiag(u(:)))' + par.beta*dmk2(:) + ...
                get_drNduT4(M,nt,dt,par,Rho(:,end) - par.drhoN) + par.gamma*dt*par.hd*get_dRudu(u,nt,par)';

    fprintf('%3d.%d\t      %3.2e \t     ||g|| = %3.2e\n',i,0,phi,norm(g));

    H13       = par.beta*2*dt*par.hd*sdiag(Rho(:)'*Abig) + par.gamma*dt*par.hd.*kron(speye(nt*par.dim),(par.Grad)'*par.Grad);
    H        = @(x) H13*x + get_drNduTdrNdu3(M,nt,dt,par,x);                           
    
    [s,pcgflag,relres,iter]    = pcg(H,-g,0.01,par.niter_pcg);
    
    dir_deriv = s'*g; % here s is the search direction, g is the gradient direction, need descent direction
    if dir_deriv >= 0
        fprintf('Warning: Search direction is not a descent direction (s''*g = %e). Reverting to gradient descent.\n', dir_deriv);
        s = -g; % back to gradient descent
        dir_deriv = s'*g;
    end

    if pcgflag ~= 0
      warning('MATLAB:pcgExitFlag','Warning: GNblock_u.m >>> iter %d, while finding s, pcg exit flag = %d \nrelres = %3.2e, iter = %d, %s',i,pcgflag,relres,iter,tag_str)
    end
    
    muls     = 0.5; % initial step size for line search
    lsiter = 1;
    while 1
        ut   = u(:) + muls*s;
        
        [phit, mk_t, phiN_t, Ru_t] = get_phi(Rho_i,ut,par);
        
        fprintf('%3d.%d\t      %3.2e \t     phit  = %3.2e        %s\n',i,lsiter,phi,phit,tag_str);
        
        % test for line search termination
        if phit < phi + 1e-4*muls*dir_deriv %1e-4*s'*g , Use Wolfe/Armijo condition with smaller muls, then this update is success
            break;                      %breaks while loop entirely (and goes to next statement after end of while loop)
        end
        muls = muls/2; lsiter = lsiter+1;
        
        % fail if small step still cannot improve, move toward a local minima or wrong direction
        if lsiter > 16
            fprintf('LSB\n');
            % ut = u;  Do not use an illegal update of ut.   
            flag = 1;
            break; % Use break instead of return to allow saving history
        end
    end
    
    % --- MODIFICATION: Record the loss components from the successful step ---
    loss_history.mk(end+1) = mk_t;
    loss_history.phiN(end+1) = phiN_t;
    loss_history.Ru(end+1) = Ru_t;
    loss_history.total(end+1) = phit;
    % --- END MODIFICATION ---

    if flag
        break; % Use break instead of return (must leave an LSB output)
    end
    u   = ut;
    phi = phit;
    
end

% --- MODIFICATION: Visualize loss history in subplots and save the figure ---
fig = figure('Visible', 'off', 'Name', 'Loss History');

% Subplot 1: MKdist Term
subplot(3, 1, 1);
plot(loss_history.mk, 'b-o');
% set clim the max and min of loss
ylim([min(loss_history.mk)*0.9-0.1, max(loss_history.mk)*1.1]);

title('MKdist Term (mk)');
ylabel('Loss');
grid on;
set(gca, 'YScale', 'log');

% Subplot 2: Image Mismatch Term
subplot(3, 1, 2);
plot(loss_history.phiN, 'r-o');
ylim([min(loss_history.phiN)*0.9-0.1, max(loss_history.phiN)*1.1]);

title('Image Mismatch Term (phiN)');
ylabel('Loss');
grid on;
set(gca, 'YScale', 'log');

% Subplot 3: Regularization Term
subplot(3, 1, 3);
plot(loss_history.Ru, 'g-o');
ylim([min(loss_history.Ru)*0.9-0.1, max(loss_history.Ru)*1.1]);
title('Regularization Term (Ru)');
xlabel('Optimization Step');
ylabel('Loss');
grid on;
set(gca, 'YScale', 'log');

% Add a main title
sgtitle(sprintf('Loss Components During Optimization: %s', tag_str), 'Interpreter', 'none');

% Create the output directory for loss plots if it doesn't exist
loss_plot_dir = fullfile(par.out_dir, 'loss');
if ~exist(loss_plot_dir, 'dir')
    mkdir(loss_plot_dir);
end

% Save the figure
saveas(fig, fullfile(loss_plot_dir, sprintf('loss_curve_%s.png', tag_str)));
close(fig);

% --- END MODIFICATION ---

end


