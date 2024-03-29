clear all;
close all;
totalgeneration = 10;
population = 50;
numcompete = 2;
beta = 1;   % Beta for mutation function

load('material_properties.mat');

figure(1);
% allfit = [];
sizeval = 300;
% Generate initial parents
for i = 1:population
    Parents(i) = beamdesign();
    [fitnesses] = getFitness(Parents(i));
    scatter3(fitnesses(1),fitnesses(2),fitnesses(3),sizeval,[0 0 1],'.');
    drawnow
    xlabel('Moment');
    ylabel('Weight');
    zlabel('Deflection');
    hold on
end

ParetoDesigns = beamdesign.empty;
for currentGeneration = 1:totalgeneration
    currentGeneration
    redval = currentGeneration/(totalgeneration + 1);
    blueval = 1-redval;
    rgbvals = [redval 0 blueval];
    winners = tournament(Parents,numcompete);
    children = [];
    disp('Creating children');
    for i = 1:length(winners)/2
        [child1,child2] = crossOver(winners(i),winners(i+1),MaterialProperties);
        children = [children,child1,child2];
    end
    disp('Mutating children');
    for i = 1:length(children)
        result(i) = mutate(children(i),currentGeneration,totalgeneration,beta,MaterialProperties);
    end
    
    
    %Elitism
    disp('Elitism');
    eliSet = [Parents,result];
    for i = 1:length(eliSet)
        FitnessOutputs(i) = maximin(i,eliSet);
    end
    keepSize = length(Parents);
    [B,I] = mink(FitnessOutputs,keepSize);
    for i = 1:keepSize
        Parents(i) = eliSet(I(i));
        [fitnesses] = getFitness(Parents(i));
        scatter3(fitnesses(1),fitnesses(2),fitnesses(3),sizeval,rgbvals,'.');
        drawnow
        if FitnessOutputs(I(i)) < 0
            ParetoDesigns(currentGeneration,i) = Parents(i);
        end
    end
    [value, index] = min(FitnessOutputs);
    [scoreHistory(currentGeneration), ind] = min(FitnessOutputs); % Save the maximum score from each generation (what does ind do here?)
%     FitnessOutputs(ind) = -Inf;
%     [secondFitness, ind2] = max(FitnessOutputs);
%     FitnessOutputs(ind2) = -Inf;
%     [thirdFitness, ind3] = max(FitnessOutputs);
    BestDesign = eliSet(ind);
%     secondDesign = eliSet(ind2);
%     thirdDesign = eliSet(ind3);
    designHistory(currentGeneration) = BestDesign;
end
value
index

