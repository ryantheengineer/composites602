function [mapped_child] = map_nplies_w(parent,child)
% Create mapping of parent to child (precedes uniform crossover)

% Rescale parent to new parameter (ex: nplies_w)
pvec = linspace(1,parent.nplies_w,parent.nplies_w);
pvec_scaled = rescale(pvec);

% Temporary child vector for mapping
tempvec = linspace(1,child.nplies_w,child.nplies_w);
tempvec_scaled = rescale(tempvec);

% Determine mapping
mapped_child = zeros();
for m = 1:length(tempvec_scaled)
    [~,mapped_child(m)] = min(abs(pvec_scaled - tempvec_scaled(m)));
end

end