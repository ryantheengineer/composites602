generation = 10;
population = 10;
numcompete = 2;
beta = 1;   % Beta for mutation function

load('material_properties.mat');

% Generate initial parents
for i = 1:population
    Parents(i) = beamdesign();
end

BestDesign = [];
secondDesign = [];
thirdDesign = [];
for currentGeneration = 1:generation
    currentGeneration
    winners = tournament(Parents,numcompete);
    children = [];
    for i = 1:length(winners)/2
        [child1,child2] = crossOver(winners(i),winners(i+1),MaterialProperties);
        children = [children,child1,child2];
    end
    for i = 1:length(children)
        result(i) = mutate(children(i),currentGeneration,generation,beta);
    end
    
    
    %% BELOW HERE IS STILL THE OLD CODE FROM DOMINIOPT
    %Elitism
    eliSet = [Parents,result,BestDesign,secondDesign,thirdDesign];
    for i = 1:length(eliSet)
        Fitness(i) = maximin(eliSet(i),eliSet);
    end
    keepSize = length(Parents);
    [B,I] = mink(Fitness,keepSize);
    for i = 1:keepSize
        Parents(i) = eliSet(I(i));
    end
    [value, index] = max(Fitness);
    [scoreHistory(currentGeneration), ind] = max(Fitness); % Save the maximum score from each generation (what does ind do here?)
    Fitness(ind) = -Inf;
    [secondFitness, ind2] = max(Fitness);
    Fitness(ind2) = -Inf;
    [thirdFitness, ind3] = max(Fitness);
    BestDesign = eliSet(ind);
    secondDesign = eliSet(ind2);
    thirdDesign = eliSet(ind3);
    playHistory(currentGeneration) = BestDesign;
end
value
index
