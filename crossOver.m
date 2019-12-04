function [child1, child2] = crossOver(parent1, parent2, MaterialProperties)
crossThres = 0.6; % not sure if this is necessary (used in uniform crossover below)

child1 = beamdesign();
child2 = beamdesign();
% NEED A WAY TO CONSTRUCT BEAMDESIGN OBJECT WITH CHOSEN PARAMETERS INSTEAD
% OF RANDOM VALUES


%% Blend crossover for nplies parameters
% Decimal blend values for nplies parameters
[child1_nplies_f1,child2_nplies_f1]     = blend_cross(parent1.nplies_f1,parent2.nplies_f1);
[child1_nplies_f2,child2_nplies_f2]     = blend_cross(parent1.nplies_f2,parent2.nplies_f2);
[child1_nplies_w,child2_nplies_w]       = blend_cross(parent1.nplies_w,parent2.nplies_w);
[child1_nplies_same,child2_nplies_same] = blend_cross(parent1.nplies_same,parent2.nplies_same);

% Round nplies parameters to the nearest even integer, minimum 2
child1.nplies_f1    = roundeven(child1_nplies_f1);
child2.nplies_f1    = roundeven(child2_nplies_f1);
child1.nplies_f2    = roundeven(child1_nplies_f2);
child2.nplies_f2    = roundeven(child2_nplies_f2);
child1.nplies_w     = roundeven(child1_nplies_w);
child2.nplies_w     = roundeven(child2_nplies_w);

% Additional steps for nplies_same parameters to make sure it's not too
% large for the given cross-section
child1_nplies_same = roundeven(child1_nplies_same);
child2_nplies_same = roundeven(child2_nplies_same);

child1_minsame = min([child1.nplies_f1,child1.nplies_f2,child1.nplies_w]);
child2_minsame = min([child2.nplies_f1,child2.nplies_f2,child2.nplies_w]);

if child1_nplies_same > child1_minsame
    child1.nplies_same = child1_minsame;
else
    child1.nplies_same = child1_nplies_same;
end

if child2_nplies_same > child2_minsame
    child2.nplies_same = child2_minsame;
else
    child2.nplies_same = child2_nplies_same;
end


%% Blend crossover for member length parameters
% Decimal blend values for nplies parameters
[child1.b_f1,child2.b_f1]   = blend_cross(parent1.b_f1,parent2.b_f1);
[child1.b_f2,child2.b_f2]   = blend_cross(parent1.b_f2,parent2.b_f2);
[child1.h_w,child2.h_w]     = blend_cross(parent1.h_w,parent2.h_w);


%% Uniform crossover for layup parameters
% Function for mapping and crossing layup
% 1. Create mapping of parent1 to child
% 2. Create mapping of parent2 to child
% 3. Uniform crossover between mapped parent1 and parent2 to child

% Do the above function for child1 and child2


% Set thickness parameters for children
get_thicknesses(child1,MaterialProperties);
get_thicknesses(child2,MaterialProperties);

% Check that member length parameters are above minimums defined by member
% thickness parameters, and if not reset them so they meet those minimums



%% Derive remaining parameter values and set them in children



%% DOMINIOPT CODE FROM HERE DOWN
% uniform cross over for all binary representation
numCards = length(parent1.gain_priority);
if (rand < crossThres)
    for i = 1:numCards
        %second row of gain priority
        if (rand <= 0.5)
            temp = parent1.gain_priority(2,i);
            parent1.gain_priority(2,i) = parent2.gain_priority(2,i);
            parent2.gain_priority(2,i) = temp;
        end
    end
end
if (rand < crossThres)
    for i = 1:numCards
        % first row of gain_cutoffs
        if (rand <= 0.5)
            temp = parent1.gain_cutoffs(1,i);
            parent1.gain_cutoffs(1,i) = parent2.gain_cutoffs(1,i);
            parent2.gain_cutoffs(1,i) = temp;
        end
    end
end
if (rand < crossThres)
    for i = 1:numCards
        %second row of trash_priority
        if (rand <= 0.5)
            temp = parent1.trash_priority(2,i);
            parent1.trash_priority(2,i) = parent2.trash_priority(2,i);
            parent2.trash_priority(2,i) = temp;
        end
    end
end
% blend cross over at each column for other representation
for i = 1:numCards
    %second row of gain cutoffs
    if (rand < crossThres)
        r = rand;
        x1 = parent1.gain_cutoffs(2,i);
        x2 = parent2.gain_cutoffs(2,i);
        parent1.gain_cutoffs(2,i) = r*x1+(1-r)*x2;
        parent2.gain_cutoffs(2,i) = (1-r)*x1+r*x2;
    end
    % third row of gain cutoffs
    if (rand < crossThres)
        r = rand;
        x1 = parent1.gain_cutoffs(3,i);
        x2 = parent2.gain_cutoffs(3,i);
        parent1.gain_cutoffs(3,i) = r*x1+(1-r)*x2;
        parent2.gain_cutoffs(3,i) = (1-r)*x1+r*x2;
    end
end
child1 = parent1;
child2 = parent2;


%% Local functions %%
function [outval] = roundeven(inval)
% Round inval to the nearest even integer
    if inval < 2
        error('roundeven function only takes values above 2');
    end
    
    extra = mod(inval,2);
    
    if extra > 1
        outval = inval + (1-(extra-1));
    else
        outval = inval - extra;
    end

end


function [child1_param,child2_param] = blend_cross(parent1_param,parent2_param)
% General blend crossover function
    r = unifrnd(0,1,1,1); % Random number to determine blend values

    child1_param = r*parent1_param + (1-r)*parent2_param;
    child2_param = (1-r)*parent1_param + r*parent2_param;

end


function get_thicknesses(child,MaterialProperties)
% Add up the total thickness of each member
    child.t_f1 = 0;
    child.t_f2 = 0;
    child.t_w = 0;

    for m = 1:child.nplies_f1
        rowname = MaterialProperties.Properties.RowNames(child.layup_f1(m,1));
        child.t_f1 = child.t_f1 + MaterialProperties.t(rowname);                
    end

    for m = 1:child.nplies_f2
        rowname = MaterialProperties.Properties.RowNames(child.layup_f2(m,1));
        child.t_f2 = child.t_f2 + MaterialProperties.t(rowname);                
    end

    for m = 1:child.nplies_w
        rowname = MaterialProperties.Properties.RowNames(child.layup_w(m,1));
        child.t_w = child.t_w + MaterialProperties.t(rowname);                
    end
end


function [child_layup] = layup_crossover(parent1,parent2,child,MaterialProperties)
% Function for mapping and crossing layup parameters (uniform crossover)

%% 1. Create mapping of parent1 to child (this might be a function in itself)
% nplies_f1
% Rescale parent to new nplies
p1_nplies_f1 = parent1.nplies_f1;
p1_nplies_f1_scaled = rescale(p1_nplies_f1);

c_nplies_f1 = child.nplies_f1;
c_nplies_f1_scaled = rescale(c_nplies_f1);

for m = 1:length(p1_nplies_f1_scaled)
    [~,ix] = min(abs(c_nplies_f1_scaled - p1_nplies_f1_scaled(m)));
end
% ix gives the index that should be used from p1_nplies_f1_scaled to map it
% to the child (assuming I wrote that correctly)

% nplies_f2

% nplies_w

% 2. Create mapping of parent2 to child


% 3. Uniform crossover between mapped parent1 and parent2 to child


% Match up columns 3 and 4 of the layup with their chosen materials



% Apply nplies_same rule

end


end