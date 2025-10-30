function [t,x,u,ux,ut,tb,S] = LR_sum_solver(Q,c,cp,G,L0,R0,tfinal,tswitch,Nt,Nx,xbuf)

% equally spaced time steps!
t = linspace(0,tfinal,Nt)';
[~,i_split] = min(abs(t-tswitch));
tmol = t(1:i_split);
ttrans = t(i_split:end);
dt = t(2)-t(1);

% Solve the potential system for r = R_x and l = L_x
par.PDE = @(h,t,u) potential_system(c,G,h,t,u);

par.ICs = {matlabFunction(diff(L0,argnames(L0))), matlabFunction(diff(R0,argnames(R0)))};
par.tspan  = tmol;
par.N      = Nx;
par.xbuf   = xbuf;
par.domain = [-par.tspan(end)-par.xbuf, par.tspan(end)+par.xbuf, par.N];
% the system of ODEs shouldn't be stiff 
par.solver = @ode78; % ode78 appears to be the fastest for this problem
par.options = odeset(RelTol=1e-10, AbsTol=1e-10);

tic
[tmol,xmol,U] = MoL1D(par.solver, par.PDE, par.tspan, par.domain, par.ICs, par.options);
Lx1 = U(:,:,1); Rx1 = U(:,:,2);
disp("Finished computing with MoL")
toc
dx = xmol(2)-xmol(1);
L1 = cumtrapz(xmol,Lx1,2);
R1 = cumtrapz(xmol,Rx1,2);
ux1 = Lx1+Rx1;
u1  = L1+R1;
ut1 = gradient(u1,dt);


% transport_solver requires the initial conditions be functions (not lists
% of discrete points), so using the solution from MoL is a bit tricky. We
% will build Chebfuns from the equally spaced data and differentiate them

% TODO try pchip, makima, or spline to see if the random spikes disappear...
% interp_style = 'linear';
interp_style = 'pchip';
Lxinterp = griddedInterpolant(xmol,Lx1(end,:),interp_style,'nearest');
Lxxinterp = griddedInterpolant(xmol,gradient(Lx1(end,:),dx),interp_style,'nearest');

Rxinterp = griddedInterpolant(xmol,Rx1(end,:),interp_style,'nearest');
Rxxinterp = griddedInterpolant(xmol,gradient(Rx1(end,:),dx),interp_style,'nearest');

% Restrict the domains of L and R
Ldomain = xmol; Rdomain = xmol;

tic
% Solve for L
[Xl,tl,l2,sl,tb_l] = transport_solver(@(ux)-Q(ux),@(ux)-c(ux),@(ux)-cp(ux), ...
    @(x)Lxinterp(x),@(x)Lxxinterp(x),Ldomain,ttrans,false);
for b = tb_l'
    tbi = find(abs(tl-b) < 1e-14);
    if(~isempty(tbi)), Xl(tbi,:) = []; tl(tbi) = []; l2(tbi,:) = []; sl(tbi,:) = []; end
end
L2 = NaN*ones(size(l2)); idxs = ~isnan(Xl(1,:));
for i = 1:length(tl)
    L2(i,idxs) = cumtrapz(Xl(i,idxs),l2(i,idxs));
end
Lt2(:,idxs) = gradient(L2(:,idxs)', dt)' - gradient(Xl(:,idxs)', dt)'.*l2(:,idxs);

% Repeat for R
[Xr,tr,r2,sr,tb_r] = transport_solver(Q,c,cp, ...
    @(x) Rxinterp(x),@(x)Rxxinterp(x),Rdomain,ttrans,false);
for b = tb_r'
    tbi = find(abs(tr-b) < 1e-14);
    if(~isempty(tbi)), Xr(tbi,:) = []; tr(tbi) = []; r2(tbi,:) = []; sr(tbi,:) = []; end
end
R2 = NaN*ones(size(r2)); idxs = ~isnan(Xr(1,:));
for i = 1:length(tr)
    R2(i,idxs) = cumtrapz(Xr(i,idxs),r2(i,idxs));
end
Rt2(:,idxs) = gradient(R2(:,idxs)', dt)' - gradient(Xr(:,idxs)', dt)'.*r2(:,idxs);
disp("Found long-time dynamics with transport_solver")
toc

tb = uniquetol([tb_l, tb_r],4*eps);

% Merge L and R. TODO This could be done far better within
% transport_solver, if it had the functionality.
x2 = cell(length(ttrans)-1,1); u2 = cell(length(ttrans)-1,1);
ux2 = cell(length(ttrans)-1,1); ut2 = cell(length(ttrans)-1,1);
for i = 1:length(ttrans)-1
    % Splitting the domain about 0 like this also removes NaNs
    idl = Xl(i+1,:) < 0;
    idr = Xr(i+1,:) > 0;
    x2{i} = [Xl(i+1,idl),Xr(i+1,idr)];
    u2{i} = [L2(i+1,idl),R2(i+1,idr)];
    ux2{i} = [l2(i+1,idl),r2(i+1,idr)];
    ut2{i} = [Lt2(i+1,idl),Rt2(i+1,idr)];
end

x = [repmat(mat2cell(xmol,1,length(xmol)),[length(tmol),1]); x2];
u = [mat2cell(u1,ones(1,length(tmol)),length(xmol)); u2];
ux = [mat2cell(ux1,ones(1,length(tmol)),length(xmol)); ux2];
ut = [mat2cell(ut1,ones(1,length(tmol)),length(xmol)); ut2];
S = [sl,sr];

end



% Ensure for u(x,t) that u(0,t) = u(L,t).
function [xfor,xint,xbac] = periodic_BCs(Nx)
    xfor = [2:Nx,1]; xint = 1:Nx; xbac = [Nx,1:Nx-1];
end

% Solve the system of PDE
%   l_t = F(l)l_x + G(l,l_x,r,r_x)
%   r_t = - F(r)r_x - G(l,l_x,r,r_x)
function u = potential_system(F,G,h,~,u)
% construct indices that ensure periodic boundary conditions are satisfied
u = reshape(u,[],2); Nx = size(u,1);
[xfor,xint,xbac] = periodic_BCs(Nx);

% Update the PDE on the interior of the domain
l = u(xint,1); lx = (u(xfor,1) - u(xbac,1))/(2*h);
r = u(xint,2); rx = (u(xfor,2) - u(xbac,2))/(2*h);

interaction = G(l,lx,r,rx);
u(xint,1) = F(l).*lx + interaction;
u(xint,2) = - F(r).*rx - interaction;

% reshape u_ut to build a long column vector as ode23 expects
u = reshape(u,Nx*2,1);
end