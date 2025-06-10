function [tsub,I] = subvector(t,delta_t)
% Return a subvector with a (larger) stepsize \approx delta_t
    tsub = linspace(t(1),t(end),(t(end)-t(1))/delta_t);
    if length(tsub) > length(t)
        tsub = t;
    end

    I = ones(size(tsub));
    for k = 1:length(tsub)
        [~,I(k)] = min(abs(t - tsub(k)));
    end
    tsub = t(I);
end