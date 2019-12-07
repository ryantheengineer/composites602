totalgeneration = 50;
population = 20;
numcompete = 3;
beta = 1;   % Beta for mutation function

load('material_properties.mat');

% Generate initial parents
for i = 1:population
    Parents(i) = beamdesign();
end

BestDesign = [];
secondDesign = [];
thirdDesign = [];
for currentGeneration = 1:totalgeneration
    currentGeneration
    winners = tournament(Parents,numcompete);
    children = [];
    for i = 1:length(winners)/2
        [child1,child2] = crossOver(winners(i),winners(i+1),MaterialProperties);
        children = [children,child1,child2];
    end
    for i = 1:length(children)
        result(i) = mutate(children(i),currentGeneration,totalgeneration,beta,MaterialProperties);
    end
    
    
    %Elitism
    eliSet = [Parents,result];
    for i = 1:length(eliSet)
        FitnessOutputs(i) = maximin(i,eliSet);
    end
    keepSize = length(Parents);
    [B,I] = mink(FitnessOutputs,keepSize);
    for i = 1:keepSize
        Parents(i) = eliSet(I(i));
    end
    [value, index] = max(FitnessOutputs);
    [scoreHistory(currentGeneration), ind] = min(FitnessOutputs); % Save the maximum score from each generation (what does ind do here?)
%     FitnessOutputs(ind) = -Inf;
%     [secondFitness, ind2] = max(FitnessOutputs);
%     FitnessOutputs(ind2) = -Inf;
%     [thirdFitness, ind3] = max(FitnessOutputs);
    BestDesign = eliSet(ind);
%     secondDesign = eliSet(ind2);
%     thirdDesign = eliSet(ind3);
    playHistory(currentGeneration) = BestDesign;
end
value
index
