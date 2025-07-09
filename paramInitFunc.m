function par_out = paramInitFunc(cfg)
% this function returns the par structure
n = cfg.true_size';
nt = cfg.nt;
dt = cfg.dt;
sigma = cfg.sigma;
add_source = cfg.add_source;
gamma = cfg.gamma;
beta = cfg.beta;
niter_pcg = cfg.niter_pcg;
dTri = cfg.dTri;
out_dir = cfg.out_dir;
maxUiter = cfg.maxUiter;
par.dim   = length(n);

% Domain size, uniform grid
n1 = n(1);
h1 = ones(n1,1);

if par.dim==2
    n2 = n(2); 
    h2 = ones(n2,1);
    [Xc, Yc] = getCellCenteredGrid(h1,h2);
    BC = {'ccn' 'ccn'};
    Grad = getCellCenteredGradMatrix(BC,h1,h2);                
    Mdis = -sigma*(Grad')*Grad;
    I = speye(prod(n));
    par.h1 = h1; par.h2 = h2;
    par.hd = h1(1)*h2(1);
    par.Xc = Xc; par.Yc = Yc;
    par.bc = getBCString(dTri);
elseif par.dim==3
    n2 = n(2); n3 = n(3);
    h2 = ones(n2,1); h3 = ones(n3,1);
    [Xc, Yc, Zc] = getCellCenteredGrid(h1,h2,h3);
    BC = {'ccn' 'ccn' 'ccn'};
    Grad = getCellCenteredGradMatrix(BC,h1,h2,h3);                
    Mdis = -sigma*(Grad')*Grad;
    I = speye(prod(n));
    par.h1 = h1; par.h2 = h2; par.h3 = h3;
    par.hd = h1(1)*h2(1)*h3(1);
    par.Xc = Xc; par.Yc = Yc; par.Zc = Zc;
    par.bc = getBCString(dTri);
else
    warning('In paramInitFunc.m: dimension of data should be either 2 or 3')
    par_out = [];
    return;
end

% Common assignments
par.add_source = add_source;
par.n     = n;
par.dt    = dt;
par.nt    = nt;
par.sigma = sigma;
par.Grad  = Grad;
par.B     = I - dt*Mdis;


if add_source
    % Placeholder for source initialization if needed
    % par.qex = ...
end

% Inversion parameters
par.gamma     = gamma;
par.alpha     = 0.1;
par.beta      = beta;
par.minUpdate = 1;
par.maxCGiter = 10;
par.niter_pcg = niter_pcg;

par.maxUiter  = maxUiter; % WARNING: more step of optimization.
par.out_dir = out_dir;

par_out = par;

end

function bc = getBCString(dTri)
    switch dTri
        case 1
            bc = 'closed';
        case 3
            bc = 'open';
        otherwise
            bc = 'closed';
    end
end


