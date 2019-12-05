function [score] = maximin(design,designpop)
% Take in a beamdesign and calculate the maximin score against the
% population of beamdesigns
population = length(designpop);
numobj = length(design.fitnesses);

for i = 1:population
    if designpop(i).fitnesses == design.fitnesses
        idx = i;
    end
end

designpop(idx) = [];

minvals = zeros((population-1),numobj);

for i = 1:population-1 % Each row is a different design
    for j = 1:numobj % Each column is a different objective
        minvals(i,j) = design.fitnesses(j) - designpop(i).fitnesses(j);   % REMEMBER TO SCALE OBJECTIVE FUNCTION VALUES TO BE ON THE SAME ORDER OF MAGNITUDE!!!!!!!!     
    end
end

mins = zeros(population-1,1);
for i = 1:length(mins)
    mins(i) = min(minvals(i,:));
end

score = max(mins);

end