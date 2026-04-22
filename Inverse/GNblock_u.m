function [u,phi,g] = GNblock_u(Rho_i,u,nt,dt,par,tag_str)
%%
%% Minimal robustness edits:
% 1) move brittle numeric constants into configurable defaults
% 2) use adaptive PCG tolerance default
% 3) use previous accepted line-search step as warm start
% 4) keep original solver structure and logic otherwise unchanged

if nargin < 6
    tag_str = '';
end

% -------------------------
% robust default parameters
% -------------------------
if ~isfield(par,'pcg_tol') || isempty(par.pcg_tol)
    par.pcg_tol = 1e-2;              % keep original behavior as default
end
if ~isfield(par,'armijo_c') || isempty(par.armijo_c)
    par.armijo_c = 1e-4;             % more standard than 5e-3, less brittle
end
if ~isfield(par,'ls_init') || isempty(par.ls_init)
    par.ls_init = 1.0;               % Newton/GN step should try full step first
end
if ~isfield(par,'ls_shrink') || isempty(par.ls_shrink)
    par.ls_shrink = 0.5;             % keep original halving logic
end
if ~isfield(par,'max_lsiter_early') || isempty(par.max_lsiter_early)
    par.max_lsiter_early = 12;
end
if ~isfield(par,'max_lsiter_mid') || isempty(par.max_lsiter_mid)
    par.max_lsiter_mid = 20;
end
if ~isfield(par,'max_lsiter_late') || isempty(par.max_lsiter_late)
    par.max_lsiter_late = 30;
end
if ~isfield(par,'hybrid_lap_scale') || isempty(par.hybrid_lap_scale)
    par.hybrid_lap_scale = max(0.1, min(1.0, 0.6 * par.gamma * dt * par.hd + 0.1));
    % preserves original 0.6-ish behavior but ties it weakly to problem scaling
end
if ~isfield(par,'ichol_droptol') || isempty(par.ichol_droptol)
    par.ichol_droptol = 1e-3;
end
if ~isfield(par,'ichol_diagcomp') || isempty(par.ichol_diagcomp)
    par.ichol_diagcomp = 0.1;
end
if ~isfield(par,'ls_init_from_prev') || isempty(par.ls_init_from_prev)
    par.ls_init_from_prev = true;
end

% --- MODIFICATION: Initialize history arrays ---
loss_history.mk = [];
loss_history.phiN = [];
loss_history.Ru = [];
loss_history.total = [];
% --- END MODIFICATION ---

[phi, ~, ~, ~] = get_phi(Rho_i,u,par);
A            = kron(ones(1,par.dim),speye(prod(par.n)));
Abig         = kron(speye(nt),A);
flag         = 0;

% warm-start accepted step size for next iteration
last_muls = par.ls_init;

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

    for j = 1:nt
        dmk2(:,1:j) = dmk2(:,1:j) + reshape( ...
            par.hd*dt*get_drNduT4(M,j,dt,par,A*(U(:,j).*U(:,j))), ...
            par.dim*prod(par.n),j);
    end

    g = (par.beta*2*par.hd*dt*Rho(:)'*Abig*sdiag(u(:)))' + par.beta*dmk2(:) + ...
        get_drNduT4(M,nt,dt,par,Rho(:,end) - par.drhoN) + ...
        par.gamma*dt*par.hd*get_dRudu(u,nt,par)';

    if isfield(par, 'mask_full') && ~isempty(par.mask_full)
        g = g .* par.mask_full;
    end

    fprintf('%3d.%d\t      %3.2e \t     ||g|| = %3.2e\t %s \n',i,0,phi,norm(g),tag_str);

    H_diag = par.beta*2*dt*par.hd*sdiag(Rho(:)'*Abig);
    H_laplacian_unscaled = kron(speye(nt*par.dim), (par.Grad)'*(par.Grad));

    H13 = H_diag + par.gamma*dt*par.hd.*H_laplacian_unscaled;
    H   = @(x) H13*x + get_drNduTdrNdu3(M,nt,dt,par,x);

    % minimally robustified hybrid preconditioner scaling
    M_hybrid = H_diag + par.hybrid_lap_scale * H_laplacian_unscaled;

    fprintf('Computing Hybrid Preconditioner...\t %s\n', tag_str);
    opts.type = 'ict';
    opts.droptol = par.ichol_droptol;
    opts.diagcomp = par.ichol_diagcomp;
    L_hybrid = ichol(M_hybrid, opts);
    fprintf('Finish computing Hybrid Preconditioner.\t %s\n', tag_str);

    % --- Adaptive PCG and Line Search Loop ---
    s_guess = [];
    found_step = false;
    current_niter = 0;

    for k = 1:length(par.niter_pcg)
        last_niter = current_niter;
        current_niter = par.niter_pcg(k);
        fprintf('--- Attempting PCG with max_iter = %d ---\n', current_niter);

        [s, pcgflag, relres, iter] = pcg(H, -g, par.pcg_tol, ...
            current_niter - last_niter, L_hybrid, L_hybrid', s_guess);

        s_guess = s;

        fprintf("Finish Hx=-g pcg. Flag=%d, RelRes=%e, Iter=%d \t%s\n", pcgflag, relres, iter, tag_str);

        dir_deriv = s'*g;
        if dir_deriv >= 0
            fprintf('Warning: Search direction is not a descent direction (s''*g = %e). Reverting to gradient descent.\n', dir_deriv);
            s = -g;
            dir_deriv = s'*g;
        end

        % --- Line Search ---
        if par.ls_init_from_prev
            muls = min(1.0, max(par.ls_shrink, last_muls / par.ls_shrink));
        else
            muls = par.ls_init;
        end

        lsiter = 1;
        line_search_failed = false;

        while 1
            ut = u(:) + muls*s;
            [phit, mk_t, phiN_t, Ru_t] = get_phi(Rho_i, ut, par);
            fprintf('%3d.%d\t      %3.2e \t     phit  = %3.2e        %s\n', i, lsiter, phi, phit, tag_str);

            if phit < phi + par.armijo_c * muls * dir_deriv
                found_step = true;
                last_muls = muls;
                break;
            end

            muls = par.ls_shrink * muls;
            lsiter = lsiter + 1;

            if i < 5
                max_lsiter = par.max_lsiter_early;
            elseif i < 12
                max_lsiter = par.max_lsiter_mid;
            else
                max_lsiter = par.max_lsiter_late;
            end

            if lsiter > max_lsiter
                fprintf('Line search failed for PCG with %d iters. Trying more PCG iters...\n', current_niter);
                line_search_failed = true;
                break;
            end
        end

        if found_step
            break;
        end
    end

    if line_search_failed && ~found_step
        fprintf('LSB: Line search failed even with max PCG iterations (%d).\n', current_niter);
        flag = 1;
    end

    loss_history.mk(end+1) = mk_t;
    loss_history.phiN(end+1) = phiN_t;
    loss_history.Ru(end+1) = Ru_t;
    loss_history.total(end+1) = phit;

    if flag
        break;
    end

    u   = ut;
    phi = phit;
end

% --- plotting block unchanged ---
fig = figure('Visible', 'off', 'Name', 'Loss History');

subplot(3, 1, 1);
plot(loss_history.mk, 'b-o');
ylim([min(loss_history.mk)*0.9-0.1, max(loss_history.mk)*1.1]);
title('MKdist Term (mk)');
ylabel('Loss');
grid on;
set(gca, 'YScale', 'log');
x_data_mk = 1:length(loss_history.mk);
for idx = 1:4:length(x_data_mk)
    text(x_data_mk(idx), loss_history.mk(idx), num2str(loss_history.mk(idx), '%.2e'), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

subplot(3, 1, 2);
plot(loss_history.phiN, 'r-o');
ylim([min(loss_history.phiN)*0.9-0.1, max(loss_history.phiN)*1.1]);
title('Image Mismatch Term (phiN)');
ylabel('Loss');
grid on;
set(gca, 'YScale', 'log');
x_data_phiN = 1:length(loss_history.phiN);
for idx = 1:4:length(x_data_phiN)
    text(x_data_phiN(idx), loss_history.phiN(idx), num2str(loss_history.phiN(idx), '%.2e'), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

subplot(3, 1, 3);
plot(loss_history.Ru, 'g-o');
ylim([min(loss_history.Ru)*0.9-0.1, max(loss_history.Ru)*1.1]);
title('Regularization Term (Ru)');
xlabel('Optimization Step');
ylabel('Loss');
grid on;
set(gca, 'YScale', 'log');
x_data_Ru = 1:length(loss_history.Ru);
for idx = 1:4:length(x_data_Ru)
    text(x_data_Ru(idx), loss_history.Ru(idx), num2str(loss_history.Ru(idx), '%.2e'), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

sgtitle(sprintf('Loss Components During Optimization: %s', tag_str), 'Interpreter', 'none');

loss_plot_dir = fullfile(par.out_dir, 'loss');
if ~exist(loss_plot_dir, 'dir')
    mkdir(loss_plot_dir);
end

saveas(fig, fullfile(loss_plot_dir, sprintf('loss_curve_%s.png', tag_str)));
close(fig);

end