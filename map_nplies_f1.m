function [mapped_child] = map_nplies_f1(parent,child)
% Create mapping of parent to child (precedes uniform crossover)

% Rescale parent to new parameter (ex: nplies_f1)
pvec = linspace(1,parent.nplies_f1,parent.nplies_f1);
pvec_scaled = rescale(pvec);

% Temporary child vector for mapping
tempvec = linspace(1,child.nplies_f1,child.nplies_f1);
tempvec_scaled = rescale(tempvec);

% Determine mapping
mapped_child = zeros();
for m = 1:length(tempvec_scaled)
    [~,mapped_child(m)] = min(abs(pvec_scaled - tempvec_scaled(m)));
end

end