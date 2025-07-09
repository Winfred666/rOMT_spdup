%function[phi,mk,phi0,phiN] =  get_phi(rho0,u,nt,dt,par)
function[phi, loss_mk, loss_phiN, loss_Ru, rho] =  get_phi(rho0, u, par)

nt = par.nt;
dt = par.dt;
rho   = advecDiff(rho0,u,nt,dt,par);

mk    = MKdist(u,rho,nt,dt, par); %=hd*dt*rho'*||v||^2
loss_mk = par.beta*mk;

loss_phiN  = 0.5*norm(rho(:,end) - par.drhoN)^2;

%% smoothing deformation field
%
Ru = 0;
if par.gamma~=0
    uvec=vec2mat(u(:),par.dim*nt);
    %GTG=par.Grad'*par.Grad;
    for ii = 1:par.dim*nt
        %Ru=Ru+0.5*par.hd*dt*(uvec(:,ii)'*GTG*uvec(:,ii));
        Ru=Ru+0.5*par.hd*dt*(norm(par.Grad*uvec(:,ii))^2); 
    end
    %}
end
loss_Ru = par.gamma.*Ru;

phi   = loss_mk + loss_phiN + loss_Ru;
end


