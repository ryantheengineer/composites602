function [child] = mutate(parentGene, currentGeneration, totalGen, beta)
numCards = length(parentGene.gain_priority);
switchThres = 0.4; % Can change this
mutateThres = 0.2; % Can change this
alpha = getalpha(currentGeneration,totalGen,beta);

allowed_angles = [0 15 30 45 60 75 90 -15 -30 -45 -60 -75]';




%% OLD DOMINIOPT CODE BELOW HERE
% for gain_priority
% for i = 1:numCards
%     % for first row (1~17, no repeated values)
%     if (rand < switchThres)
%         num1 = randi([1 numCards]);
%         num2 = parentGene.gain_priority(1,i);
%         index = find(parentGene.gain_priority(1,:)==num1, 1, 'first');
%         parentGene.gain_priority(1,index) = num2;
%         parentGene.gain_priority(1,i) = num1;      
%     end
%     % for second row mutation(binary)
%     if (rand < mutateThres)
%         parentGene.gain_priority(2,i) = round(rand);
%     end   
% end
% for gain_cutoffs
for i = 1:numCards
    % for first row (binary)
    if (rand < mutateThres)
        parentGene.gain_cutoffs(1,i) = round(rand);
    end
    % for second row (continusous number between 0 and 1)
    if (rand < mutateThres)
        r = rand;
        x = parentGene.gain_cutoffs(2,i);
        if (r <= x)
            parentGene.gain_cutoffs(2,i) = (r^alpha)*x^(1-alpha);
        elseif (r > x)
            parentGene.gain_cutoffs(2,i) = 1-((1-r)^alpha)*(1-x)^(1-alpha);
        end
%%%%%%%%%%%%%% TESTING CODE remove afterward %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if (parentGene.gain_cutoffs(2,i) > 1 || parentGene.gain_cutoffs(2,i) < 0)
%             disp(parentGene.gain_cutoffs(2,i))
%         end
            
    end
    % for third row (each col has desire max values) => later  
        
end
% for play_priority  
for i = 1:10
    if (rand < switchThres)
        num1 = randi([1 10]);
        num2 = parentGene.play_priority(1,i);
        index = find(parentGene.play_priority(1,:)==num1, 1, 'first');
        parentGene.play_priority(1,index) = num2;
        parentGene.play_priority(1,i) = num1;      
    end    
end
% for trash_priority
for i = 1:numCards
    if (rand < switchThres)
        num1 = randi([1 numCards]);
        num2 = parentGene.trash_priority(1,i);
        index = find(parentGene.trash_priority(1,:)==num1, 1, 'first');
        parentGene.trash_priority(1,index) = num2;
        parentGene.trash_priority(1,i) = num1;
    end
end

%%% gain_cutoffs thrid row
for i = 1:3
    if (rand < mutateThres)
        r = randi([2 12]);
        x = parentGene.gain_cutoffs(3,i);
        if (r <= x)
            parentGene.gain_cutoffs(3,i) = round((r^alpha)*x^(1-alpha));
        elseif (r > x)
            parentGene.gain_cutoffs(3,i) = round(12-((12-r)^alpha)*(12-x)^(1-alpha));
        end
    end
end
for i = 4:6
    if (rand < mutateThres)
        r = randi([2 12]);
        x = parentGene.gain_cutoffs(3,i);
        if (r <= x)
            parentGene.gain_cutoffs(3,i) = round((r^alpha)*x^(1-alpha));
        elseif (r > x)
            parentGene.gain_cutoffs(3,i) = round(12-((12-r)^alpha)*(12-x)^(1-alpha));
        end
    end
end
for i = 7:numCards
    if (rand < mutateThres)
        r = randi([2 12]);
        x = parentGene.gain_cutoffs(3,i);
        if (r <= x)
            parentGene.gain_cutoffs(3,i) = round((r^alpha)*x^(1-alpha));
        elseif (r > x)
            parentGene.gain_cutoffs(3,i) = round(12-((12-r)^alpha)*(12-x)^(1-alpha));
        end
    end
end
child = parentGene;

end

%% Local functions %%
function [alpha] = getalpha(currentGeneration,totalGen,beta)
% Compute alpha value for dynamic mutation
    alpha = (1-(currentGeneration-1)/totalGen)^beta;
end

