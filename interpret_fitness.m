function [FSfitnesses] = interpret_fitness(fitnesses,moment_scale,weight_scale,dmax_scale)
% Interpret the first and third fitness values as factors of safety at true
% scale

moment_objective = fitnesses(1);
weight_objective = fitnesses(2);
dmax_objective   = fitnesses(3);


moment_objective = moment_objective*moment_scale;
weight_objective = weight_objective*weight_scale;
dmax_objective   = dmax_objective*dmax_scale;

% Convert dmax to inches
dmax_objective = 0.0254*dmax_objective;
dmax_objective = 1/dmax_objective;

moment_objective = 17460/moment_objective;

FSfitnesses = [moment_objective; weight_objective; dmax_objective];


end