function [score] = maximin(designidx,designpop)
% Take in a beamdesign and calculate the maximin score against the
% population of beamdesigns
population = length(designpop);
designfit = getFitness(designpop(designidx));
numobj = size(designfit,1);

% Calculate the fitnesses and fill a 3D array of fitness values
popfitnesses = zeros(size(designfit,1),size(designfit,2),population);
for i = 1:population
    popfitnesses(:,:,i) = getFitness(designpop(i));
end

popfitnesses(:,:,designidx) = []; % Remove the design we are interested in

minvals = zeros((population-1),numobj);

for i = 1:population-1 % Each row is a different design
    for j = 1:numobj % Each column is a different objective
        minvals(i,j) = designfit(j) - popfitnesses(j,1,i);
    end
end

mins = zeros(population-1,1);
for i = 1:length(mins)
    mins(i) = min(minvals(i,:));
end

score = max(mins);

end