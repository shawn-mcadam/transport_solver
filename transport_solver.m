function [X,t,U,S,tbreaks,partition,XI] = transport_solver(Q,c,cp,phi,phip,x0,t,multivalued)
% transport_solver - Solve a 1st-order quasilinear PDE of the form
% u_t + c(u)u_x = 0 using the method of characteristics.
% 
% Syntax
%   [X,t,U,S,tbreaks,partition,XI] = transport_solver(Q,c,cp,phi,phip,x0,t,multivalued)
%
% Input Arguments
%   Q,c,cp - Flux and its derivatives
%   phi,phip - Initial condition and its derivative
%   x0 - Dictates domain boundaries & total number of space steps
%   t - Time mesh
%   multivalued - If true then solver does not check for or correct breaks
%
% Output Arguments
%   X - Nonuniform spatial mesh
%   t - Time mesh with break times appended
%   U - 2D array dictating solution values on each mesh point: U(time, space)
%   S - Locations of each break (NaN until the wave breaks): S(time,j-th break)
%   tbreaks - Times each break occurs
%   partition - S along with leftmost and rightmost characteristics
%   XI - Initial values for each characteristic curve

% TODO validate the user's arguments: Q, c, cp, phi, and phip must all be
% functions taking a single argument... multivalued is a bool.
% TODO find every sign change in c(phi(x))*phip(x) automatically (maybe
% just make each function argument into a chebfun?)
% TODO if we allow tracking special values then they should be handled
% the same way as confluence

if(~issorted(t))
    error("burgersolver given unsorted time vector t. t must be sorted")
end
t1 = t(1); t = t-t1;
if(isrow(t)), t = t'; end
if(iscolumn(x0)), x0 = x0'; end
% if(iscolumn(pts)), pts = pts'; end

if ~multivalued
    [xi,tbreaks] = burger_breakers(x0,cp,phi,phip);
    % remove any tbreaks that are after t(end)
    I = tbreaks > t(end); tbreaks(I) = []; xi(I) = [];
end

if multivalued || isempty(tbreaks)
    % We need not worry about the wave breaking, so do upwind/downwind with
    % _the_ optimal nonuniform spatial mesh

    % TODO Solve this problem the same way as the case with breaks
    % (rootfind to find the exact solution on chebyshev nodes at each time)
    % The solution is always unique because there are no breaks
    U = phi(x0);
    XI = repmat(x0,length(t),1);
    X = x0 + c(U).*t;

    % U is constant on characteristics and doesn't break:
    U = repmat(U,length(t),1); 
    partition = X(:,[1,end]); tbreaks = double.empty; S = double.empty;
    return
end


% --- Make the solution single valued & tightly bracket discontinuities ---

% append and sort tb into t and remember the position of each tb in t.
[t,I,~] = uniquetol([tbreaks;t],4*eps); k_tb(I) = 1:length(I); k_tb = k_tb(1:length(tbreaks));
[~,ia,ic] = uniquetol(tbreaks,4*eps);
k_tb = k_tb(ia(ic));


% The characteristic starting at xi(j) will end at the break s_j(tb_j).
% Find each s_j(tb_j) and construct initial guesses for following s(t)
% values.
% initialize init_guess by rootfinding equal_areas extremely carefully.
% 1. Parametrize the zero contours of line_fun by theta in [pi/2,pi]
% 2. Rootfind equal_areas along that zero contour
m = length(t); n = length(xi);
partition = zeros(m,n); init_guess = zeros(m,2,n);
for j = 1:n
    partition(1:k_tb(j),j) = xi(j) + c(phi(xi(j))).*t(1:k_tb(j));
    xi0 = [xi(j), xi(j)];
    init_guess(1:k_tb(j),:,j) = repmat(xi0,k_tb(j),1);

    % dchars_zero_contour = @(theta) xi0 + fzero(@(s) secant(t(k_tb(j)+1), xi0 + s.*[cos(theta); sin(theta)]),[0,2]) .* [cos(theta); sin(theta)];
    dchars_zero_contour = @(theta) xi0 + abs(fzero(@(s) secant(t(k_tb(j)+1), xi0 + abs(s).*[cos(theta), sin(theta)]),1e-10)) .* [cos(theta), sin(theta)];
    theta_viability = @(theta) equal_areas(dchars_zero_contour(theta));
    theta = fzero(theta_viability, [pi/2,pi]);
    xi_j = dchars_zero_contour(theta);

    init_guess(k_tb(j)+1,:,j) = xi_j;
end


fsolve_opts = optimoptions('fsolve','Display','off', ...
    FunctionTolerance=2.5e-15, ...
    OptimalityTolerance=2.5e-15, ...
    StepTolerance=2.5e-15, ...
    Algorithm='trust-region-dogleg', ...
    SpecifyObjectiveGradient=true);

for k = min(k_tb)+1:m
    % for each point of discontinuity ...
    for j = intersect(find(tbreaks' < t(k),inf), 1:n)
        xi0 = init_guess(k,:,j);
        % if any(isnan(xi0)), continue; end
        xi = fsolve(@(xi1_xi2) equal_areas_function(xi1_xi2,t(k)), xi0, fsolve_opts);
        partition(k,j) = sum(xi+c(phi(xi)).*t(k))/2;

        % the initial guess becomes initial guess bounds for this region
        init_guess(k:k+1,:,j) = [xi;xi];
    end

    % for each consecutive nonincreasing pair in the partition ...
    [~,I]=sort(partition(k,:));
    for j=(find(diff(I)<0))
        % if the break is impossible, then remove its column
        z = find(k < k_tb(j:j+1))-1; % this is an offset: either 0 or 1
        if(~isempty(z))
            partition(:,j+z) = [];
            tbreaks(j+z) = [];
            k_tb(j+z) = [];
            init_guess(:,:,j+z) = [];
            n = n-1;
        else % This is a point of confluence so init_guess must jump
            % TODO find the exact point of confluence and append that time
            xi0 = [init_guess(k,1,j),init_guess(k,2,j+1)];
            xi = fsolve(@(xi1_xi2) equal_areas_function(xi1_xi2,t(k)), xi0, fsolve_opts);
            partition(k,j:j+1) = sum(xi+c(phi(xi)).*t(k))/2;
            init_guess(k:k+1,:,j) = [xi;xi]; init_guess(k:k+1,:,j+1) = [xi;xi];

            % % Set the rightmost of the two breaks to NaN to indicate merge with a NaN
            % init_guess(k:end,:,j+1) = NaN; partition(k:end,j+1) = NaN;
            % % Move these to the end of the list
            % init_guess = init_guess(:,:,[1:j,j+2:end,j+1]); partition = partition(:,:,[1:j,j+2:end,j+1]);
            % n = n-1;
        end
    end
end
init_guess(m+1,:) = [];

% Put the list back in order and set the NaN values of partition & init_guess
% [~,I] = sort(partition(1,:));
% partition = partition(:,I);
% init_guess = init_guess(:,:,I);


% Get each break's position S at each time t from the partition
S = NaN*ones(m,n); for j=1:n, S(k_tb(j):end,j) = partition(k_tb(j):end,j); end


% figure
% plot(partition,t,Color="r"); hold on
% yline(tbreaks)
% plot(init_guess,t);
% pause(0.01)

% % plot(init_guess(:,2*j-1:2*j)+c(phi(init_guess(:,2*j-1:2*j))).*t,t)
% plot(init_guess+c(phi(init_guess)).*t,t)



% Compute u at each mesh point by tracing each (x,t) along a characteristic
% until the point (xstar,0) at which point u(x,t)=phi(xstar).
% Note that defining the endpoints of partition
% like this _guarantees_ the bracket works at the endpoints

% A root xstar of @(x) reverse(x,a,t) is the initial condition of a
% characteristic curve equal to a at time t
reverse = @(x,a,t) x + c(phi(x)).*t - a;


% TODO Do not touch X after any two breaks merge
% TODO how do I adjust the initial guess to increase a small domain's size?
% x0(1) = min(x0(1),min(partition,[],'all')-128*eps); x0(end) = max(x0(end),max(partition,[],'all')+128*eps);
partition = [x0(1) + c(phi(x0(1)))*t, partition, x0(end) + c(phi(x0(end)))*t];
init_guess = [x0(1)+0*t, init_guess, x0(end)+0*t];

NPart=ceil(length(x0)/(n+1));
X = NaN*ones(m,(n+1)*NPart); XI = X;
for j=1:n+1
    for k=1:m
        if(partition(k,j+1) - partition(k,j) < 2*eps(partition(k,j+1)))
            break
        end

        % Bounds for fzero
        bounds = [init_guess(k,2*j-1),init_guess(k,2*j)-neps(init_guess(k,2*j))];

        % Spatial mesh built from the partition
        X(k,1+(j-1)*NPart:j*NPart-1) = chebynodes(partition(k,j), partition(k,j+1), NPart-1);
        X(k,[1+(j-1)*NPart,j*NPart-1]) = bounds+c(phi(bounds))*t(k);

        % Characteristic's initial values
        XI(k,1+(j-1)*NPart+1:j*NPart-2) = arrayfun(@(a) fzero(@(xi0) reverse(xi0,a,t(k)), bounds), X(k,1+(j-1)*NPart+1:j*NPart-2) );
        XI(k,[1+(j-1)*NPart,j*NPart-1]) = bounds;
    end
end
X(:,end) = x0(end) + c(phi(x0(end)))*t;
XI(:,end) = x0(end)+0*t;


% % Sort in special values (e.g. roots and critical points)
% if ~isempty(pts)
%     X = [X,pts+c(phi(pts))*t]; XI = [XI, pts+0*t];
%     [X,I] = sort(X,2);
% end
% for k = 1:m, XI(k,:) = XI(k,I(k,:)); end


% TODO probably best to manually set U to NaN at the discontinuities
U = phi(XI);

% restore t
t = t+t1;



function x = chebynodes(a,b,n)

    x = (a+b)/2 - (b-a)*cos(pi*(2*(1:n)-1)/2/n)/2;

end

% This function's roots yield a break s(t), along with xi1, xi2 from the
% two characteristics that intersect at s(t).
% TODO implement with the symbolic toolbox's vpa and solve with quadmath
function [F,J] = equal_areas_function(xi,time)
    xi_1 = xi(1); xi_2 = xi(2);
    c1 = c(phi(xi_1)); c2 = c(phi(xi_2));
    q1 = Q(phi(xi_1)); q2 = Q(phi(xi_2));

    F = [secant(time,xi); equal_areas(xi)];

    % Compute the Jacobian
    if nargout > 1
        J = zeros(2,2);
        J(1,:) = dsecant(time,xi);
        J(2,1) = q2 - q1 + (xi_1 - xi_2)*(c1*phip(xi_1) - c1*phip(xi_1) + cp(phi(xi_1))*phi(xi_1)*phip(xi_1)) + c1*phi(xi_1) - c2*phi(xi_2) - phi(xi_1)*(c1 - c2) - cp(phi(xi_1))*integral(phi, xi_2, xi_1)*phip(xi_1);
        J(2,2) = q1 - q2 - (xi_1 - xi_2)*(c2*phip(xi_2) - c2*phip(xi_2) + cp(phi(xi_2))*phi(xi_2)*phip(xi_2)) - c1*phi(xi_1) + c2*phi(xi_2) + phi(xi_2)*(c1 - c2) + cp(phi(xi_2))*integral(phi, xi_2, xi_1)*phip(xi_2);
    end
end

function Eq = equal_areas(xi)
    xi_1 = xi(1); xi_2 = xi(2);
    c1 = c(phi(xi_1)); c2 = c(phi(xi_2));
    q1 = Q(phi(xi_1)); q2 = Q(phi(xi_2));
    Eq = (xi_1-xi_2)*(q2 - q1 - (c2*phi(xi_2) - c1*phi(xi_1))) - (c1-c2)*integral(phi,xi_2,xi_1);
end

function L = secant(time,xi)
    L = secant_fun(1e-14,@(x)c(phi(x)),@(x)cp(phi(x))*phip(x),time,xi(1),xi(2))-2*eps;
end

function dL = dsecant(time,xi)
    dL = dsecant_fun(1e-14,@(x)c(phi(x)),@(x)cp(phi(x))*phip(x),@(x)0,time,xi(1),xi(2));
end

function L = secant_fun(epsilon,f,fp,t,xi_1,xi_2)
    L = 1 + (f(xi_1)-f(xi_2))*t./(xi_1-xi_2);
    idxs    = abs(xi_1-xi_2) < epsilon;
    L(idxs) = 1 + fp(xi_1(idxs))*t;
end

function J = dsecant_fun(epsilon,f,fp,fpp,t,xi_1,xi_2)
    if abs(xi_1-xi_2) < epsilon
        J = 1 + fpp(xi_1)*t;
    else
    J(1,1) =  (fp(xi_1)*(xi_1-xi_2)-(f(xi_1)-f(xi_2)))*t./((xi_1-xi_2).^2);
    J(1,2) = -(fp(xi_2)*(xi_1-xi_2)-(f(xi_1)-f(xi_2)))*t./((xi_1-xi_2).^2);
    end
end

% taken from https://stackoverflow.com/a/50256298
function out=neps(in)
out=eps(in-eps(in));
end

end
