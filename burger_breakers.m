function [xi,tb] = burger_breakers(x0,c_prime,phi,phi_prime)
% x0 must be a vector of x0 values and at least 2 of them must bracket
% a break. Returns the earliest breaking time

% find every interval with a positive breaking time
objective = @(x) c_prime(phi(x)).*phi_prime(x);

% NOTE! it is important that xi is sorted

% TODO need another tol parameter for this... also check if the length of
% each disjoint interval is significant. 
% --- TODO maybe I can remove this tol because we remove any tbreaks greater
% than t(end)!!! -------- Can also check output before returning to see if
% tb is not too large????
i = find(objective(x0(2:end-1)) < -1e-8)+1;
% i = find(objective(x0(2:end-1)) < 0)+1;
brackets = x0(setdiff(i,intersect(i-1,i+1)));

if isempty(brackets)
    xi = NaN;
    tb = inf;
    return
end

nb = length(brackets);
if mod(nb,2)~=0
    disp("Fail in burgersolver burger_breakers. Predicting")
    disp(nb)
    disp("breaks. Pass in x0 with more points")
end

nbreaks = length(brackets)/2;
xi = zeros(1,nbreaks); tb = zeros(nbreaks,1);
for j = 1:nbreaks
    [xi(j),slope_of_caustic] = fminbnd(objective, brackets(2*j-1),brackets(2*j));
    tb(j) = - 1/slope_of_caustic;
end
[xi,I] = sort(xi); tb = tb(I);
end