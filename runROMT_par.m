function [cfg,flag] = runROMT_par(cfg)

<<<<<<< HEAD
reInitializeU = 1; %1 if reinitialize u to 0 before each time step; 0 if not, unless first time step
=======
% reInitializeU = 1; %1 if reinitialize u to 0 before each time step; 0 if not, unless first time step
>>>>>>> 9749637 (init)

if ~exist(cfg.out_dir,'dir')
    mkdir(cfg.out_dir)
end

fname = sprintf('%s/record.txt',cfg.out_dir);
if ~exist(sprintf('%s/record.txt',cfg.out_dir),'file')
    csvwrite_with_headers(fname,[0 0 0 0 0 0 0 0 0],{'time-ind','ti','tf','phi','mk','Ru','phiN','max(u)','toc'});
end

rho_n = cfg.vol(1).data(:);

<<<<<<< HEAD
if ~exist(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'file')
    save(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'rho_n');
end
=======
% if ~exist(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'file')
% Always save rho_n no matter how.
save(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'rho_n');
% end
>>>>>>> 9749637 (init)

fprintf('\n =============== rOMT Starts ===============\n')
fprintf('______________________________________________\n\n')
fprintf(' tag:\t\t%s\n dataset:\t%s\n sigma:\t\t%.4f\n gamma:\t\t%.4f\n beta:\t\t%.4f\n nt:\t\t%d\n dt:\t\t%.2f\n pcg:\t\t%d\n',cfg.tag,cfg.dataset_name,cfg.sigma,cfg.gamma,cfg.beta,cfg.nt,cfg.dt,cfg.niter_pcg)
fprintf(' size:\t\t%s\n do_resize:\t%d\n resize_factor:\t%.2f\n start frame:\t%d\n end frame:\t%d\n frame jump:\t%d\n\n\n',sprintf('%d ', cfg.true_size),cfg.do_resize,cfg.size_factor,cfg.first_time,cfg.last_time+cfg.time_jump,cfg.time_jump)
%%
%{
tag = cfg.tag;
<<<<<<< HEAD
data_dir = cfg.data_dir;
first_time = cfg.first_time;
time_jump = cfg.time_jump;
extension = cfg.extension;
=======
first_time = cfg.first_time;
time_jump = cfg.time_jump;
>>>>>>> 9749637 (init)
x_range = cfg.x_range;
y_range = cfg.y_range;
z_range = cfg.z_range;
do_resize = cfg.do_resize;
size_factor = cfg.size_factor;
smooth = cfg.smooth;
true_size = cfg.true_size;
nt = cfg.nt;
dt = cfg.dt;
sigma = cfg.sigma;
add_source = cfg.add_source;
gamma = cfg.gamma;
beta = cfg.beta;
niter_pcg = cfg.niter_pcg;
dTri = cfg.dTri;
out_dir = cfg.out_dir;
%}
%%
profile on
clear T; T = 0;
<<<<<<< HEAD
parfor (tind = 1:length(cfg.first_time:cfg.time_jump:cfg.last_time),2)
=======
% remove the worker limit to leverage more cores
global_steps = length(cfg.first_time:cfg.time_jump:cfg.last_time);
% DEBUG: change parfor to for if any problem inside.
for tind = 1:global_steps
>>>>>>> 9749637 (init)
    fprintf('tind = %d\n',tind)
    tic
    %{
    if exist(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'file')==2 && exist(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'file') == 2        
        rho_n = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind));
        rho_n = rho_n.rho_n;
        u = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind));
        u = reshape(u.u,[],cfg.nt);
        continue
    end
    %}
<<<<<<< HEAD
    rho_0 = cfg.vol(tind).data(:);

    %true final density
    par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.add_source,cfg.gamma,cfg.beta,cfg.niter_pcg,cfg.dTri);
=======
    % must load from original data.
    rho_0 = cfg.vol(tind).data(:);

    %true final density
    par = paramInitFunc(cfg);
>>>>>>> 9749637 (init)
    par.drhoN     = cfg.vol(tind+1).data(:);
    
    %par = paramInitFunc(true_size',nt,dt,sigma,add_source,gamma,beta,niter_pcg,dTri);
    %par.drhoN     = rho_N(:);
    
<<<<<<< HEAD
    % initial guess for u:
    if tind == 1 || reInitializeU
    u = zeros(par.dim*prod(par.n),par.nt);%zeros(2*prod(par.n),par.nt);
    end
    
=======
    % must reinitialize u, cannot leverage last pass:
    u = zeros(par.dim*prod(par.n),par.nt);%zeros(2*prod(par.n),par.nt);

>>>>>>> 9749637 (init)
    %% Descent for u
    fprintf('\n =============== Descent on u ===============\n')
    fprintf('______________________________________________\n\n')
    fprintf('i.lsiter\tphi    \t      descent output\n')
    fprintf('________    ___________     __________________\n')
<<<<<<< HEAD
    [u,phi,dphi] = GNblock_u(rho_0,u,par.nt,par.dt,par);
    
    
    [phi,mk,phiN,Rho,Ru]  = get_phi(rho_0,u,par.nt,par.dt,par);
    rho_n = Rho(:,end);
    btoc = toc;
    T = T + btoc;
=======
    [u,~,~] = GNblock_u(rho_0,u,par.nt,par.dt,par, sprintf("f_%d_%d",cfg.first_time + tind*cfg.time_jump, cfg.first_time + (tind+1)*cfg.time_jump));
    
    [phi,mk,phiN,Ru, Rho]  = get_phi(rho_0,u,par);
    fprintf('________    ###########     __________________\n')
    rho_n = Rho(:,end);
    btoc = toc;
    T = T + btoc;
    fprintf('########    ###########     __________________\n')
>>>>>>> 9749637 (init)
    
    dlmwrite(fname,[tind,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,phi,mk,Ru,phiN,max(u(:)),btoc],'-append');
    
    save_un(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),u);
    save_rhon(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),rho_n);
    
    fprintf('tind = %d, max(u) = %5.4f\n',tind,max(u));
end
<<<<<<< HEAD
fprintf('\n =============== rOMT Ends ===============\n')
fprintf('\n Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

profile viewer
=======

fprintf('\n =============== rOMT Ends ===============\n')
fprintf('\n Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

if usejava("desktop")
profile viewer
end
% If not using MATLAB desktop, just print the profile results to command window
>>>>>>> 9749637 (init)
profile off
flag = 1;
end