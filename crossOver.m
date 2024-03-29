function [child1, child2] = crossOver(parent1, parent2, MaterialProperties)
crossThres = 0.5; % not sure if this is necessary (used in uniform crossover below)

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

child1_maxsame = min([child1.nplies_f1,child1.nplies_f2,child1.nplies_w]);
child2_maxsame = min([child2.nplies_f1,child2.nplies_f2,child2.nplies_w]);

if child1_nplies_same > child1_maxsame
    child1.nplies_same = child1_maxsame;
else
    child1.nplies_same = child1_nplies_same;
end

if child2_nplies_same > child2_maxsame
    child2.nplies_same = child2_maxsame;
else
    child2.nplies_same = child2_nplies_same;
end

%% Uniform crossover for layup parameters
[c1_layup_f1,c1_layup_f2,c1_layup_w,c2_layup_f1,c2_layup_f2,c2_layup_w] = layup_crossover(parent1,parent2,child1,child2,crossThres);
child1.layup_f1 = c1_layup_f1;
child1.layup_f2 = c1_layup_f2;
child1.layup_w  = c1_layup_w;
child2.layup_f1 = c2_layup_f1;
child2.layup_f2 = c2_layup_f2;
child2.layup_w  = c2_layup_w;

% Set thickness parameters for children
[t_f1,t_f2,t_w] = get_thicknesses(child1,MaterialProperties);
child1.t_f1 = t_f1;
child1.t_f2 = t_f2;
child1.t_w  = t_w;
[t_f1,t_f2,t_w] = get_thicknesses(child2,MaterialProperties);
child2.t_f1 = t_f1;
child2.t_f2 = t_f2;
child2.t_w  = t_w;


%% Blend crossover for member length parameters
% Decimal blend values for nplies parameters
[child1.b_f1,child2.b_f1]   = blend_cross(parent1.b_f1,parent2.b_f1);
[child1.b_f2,child2.b_f2]   = blend_cross(parent1.b_f2,parent2.b_f2);
[child1.h_w,child2.h_w]     = blend_cross(parent1.h_w,parent2.h_w);

% Check that member length parameters are above minimums defined by member
% thickness parameters, and if not reset them so they meet those minimums
if child1.b_f1 < (child1.t_w + 0.001)
    child1.b_f1 = child1.t_w + 0.001;
end

if child1.b_f2 > (child1.t_w + 0.001)
    child1.b_f2 = child1.t_w + 0.001;
end

if child2.b_f1 < (child2.t_w + 0.001)
    child2.b_f1 = child2.t_w + 0.001;
end

if child2.b_f2 > (child2.t_w + 0.001)
    child2.b_f2 = child2.t_w + 0.001;
end


%% Derive remaining parameter values and set them in children
[A_f1,A_f2,A_w,A,ybar,d_f1,d_f2,d_w,I_f1,I_f2,I_w,I] = derive_properties(child1);
child1.A_f1 = A_f1;
child1.A_f2 = A_f2;
child1.A_w  = A_w;
child1.A    = A;
child1.ybar = ybar;
child1.d_f1 = d_f1;
child1.d_f2 = d_f2;
child1.d_w  = d_w;
child1.I_f1 = I_f1;
child1.I_f2 = I_f2;
child1.I_w  = I_w;
child1.I    = I;

[A_f1,A_f2,A_w,A,ybar,d_f1,d_f2,d_w,I_f1,I_f2,I_w,I] = derive_properties(child2);
child2.A_f1 = A_f1;
child2.A_f2 = A_f2;
child2.A_w  = A_w;
child2.A    = A;
child2.ybar = ybar;
child2.d_f1 = d_f1;
child2.d_f2 = d_f2;
child2.d_w  = d_w;
child2.I_f1 = I_f1;
child2.I_f2 = I_f2;
child2.I_w  = I_w;
child2.I    = I;


end



%%%%%%%%%%%%%%%%%%%%%
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
    r = unifrnd(0,1); % Random number to determine blend values

    child1_param = r*parent1_param + (1-r)*parent2_param;
    child2_param = (1-r)*parent1_param + r*parent2_param;

end


function [t_f1,t_f2,t_w] = get_thicknesses(child,MaterialProperties)
% Add up the total thickness of each member
%     t_f1 = 0;
%     t_f2 = 0;
%     t_w = 0;
    
    tply = MaterialProperties.t(1);
    
    t_f1 = child.nplies_f1*tply;
    t_f2 = child.nplies_f2*tply;
    t_w  = child.nplies_w*tply;
%     for m = 1:child.nplies_f1
%         rowname = MaterialProperties.Properties.RowNames(child.layup_f1(m,1));
%         t_f1 = t_f1 + MaterialProperties.t(rowname);                
%     end
% 
%     for m = 1:child.nplies_f2
%         rowname = MaterialProperties.Properties.RowNames(child.layup_f2(m,1));
%         t_f2 = t_f2 + MaterialProperties.t(rowname);                
%     end
% 
%     for m = 1:child.nplies_w
%         rowname = MaterialProperties.Properties.RowNames(child.layup_w(m,1));
%         t_w = t_w + MaterialProperties.t(rowname);                
%     end
end


function [c1_layup_f1,c1_layup_f2,c1_layup_w,c2_layup_f1,c2_layup_f2,c2_layup_w] = layup_crossover(parent1,parent2,child1,child2,crossThres)
% Function for mapping and crossing layup parameters (uniform crossover)

    %% 1. Create mapping of parent1 to children (this might be a function in itself)
    % nplies_f1
    mapped_p1c1_f1 = map_nplies_f1(parent1,child1);
    mapped_p1c2_f1 = map_nplies_f1(parent1,child2);

    % nplies_f2
    mapped_p1c1_f2 = map_nplies_f2(parent1,child1);
    mapped_p1c2_f2 = map_nplies_f2(parent1,child2);

    % nplies_w
    mapped_p1c1_w = map_nplies_w(parent1,child1);
    mapped_p1c2_w = map_nplies_w(parent1,child2);

    %% 2. Create mapping of parent2 to children
    % nplies_f1
    mapped_p2c1_f1 = map_nplies_f1(parent2,child1);
    mapped_p2c2_f1 = map_nplies_f1(parent2,child2);

    % nplies_f2
    mapped_p2c1_f2 = map_nplies_f2(parent2,child1);
    mapped_p2c2_f2 = map_nplies_f2(parent2,child2);

    % nplies_w
    mapped_p2c1_w = map_nplies_w(parent2,child1);
    mapped_p2c2_w = map_nplies_w(parent2,child2);

    %% 3. Create mapped parent and child layups
    % layup_f1
    p1c1_layup_f1 = cell(length(mapped_p1c1_f1),4);
    for j = 1:length(mapped_p1c1_f1)
        p1c1_layup_f1(j,:) = parent1.layup_f1(mapped_p1c1_f1(j),:);    
    end

    p2c1_layup_f1 = cell(length(mapped_p2c1_f1),4);
    for j = 1:length(mapped_p2c1_f1)
        p2c1_layup_f1(j,:) = parent2.layup_f1(mapped_p2c1_f1(j),:);    
    end

    p1c2_layup_f1 = cell(length(mapped_p1c2_f1),4);
    for j = 1:length(mapped_p1c2_f1)
        p1c2_layup_f1(j,:) = parent1.layup_f1(mapped_p1c2_f1(j),:);    
    end

    p2c2_layup_f1 = cell(length(mapped_p2c2_f1),4);
    for j = 1:length(mapped_p2c2_f1)
        p2c2_layup_f1(j,:) = parent2.layup_f1(mapped_p2c2_f1(j),:);    
    end


    % layup_f2
    p1c1_layup_f2 = cell(length(mapped_p1c1_f2),4);
    for j = 1:length(mapped_p1c1_f2)
        p1c1_layup_f2(j,:) = parent1.layup_f2(mapped_p1c1_f2(j),:);    
    end

    p2c1_layup_f2 = cell(length(mapped_p2c1_f2),4);
    for j = 1:length(mapped_p2c1_f2)
        p2c1_layup_f2(j,:) = parent2.layup_f2(mapped_p2c1_f2(j),:);    
    end

    p1c2_layup_f2 = cell(length(mapped_p1c2_f2),4);
    for j = 1:length(mapped_p1c2_f2)
        p1c2_layup_f2(j,:) = parent1.layup_f2(mapped_p1c2_f2(j),:);    
    end

    p2c2_layup_f2 = cell(length(mapped_p2c2_f2),4);
    for j = 1:length(mapped_p2c2_f2)
        p2c2_layup_f2(j,:) = parent2.layup_f2(mapped_p2c2_f2(j),:);    
    end


    % layup_w
    p1c1_layup_w = cell(length(mapped_p1c1_w),4);
    for j = 1:length(mapped_p1c1_w)
        p1c1_layup_w(j,:) = parent1.layup_w(mapped_p1c1_w(j),:);    
    end

    p2c1_layup_w = cell(length(mapped_p2c1_w),4);
    for j = 1:length(mapped_p2c1_w)
        p2c1_layup_w(j,:) = parent2.layup_w(mapped_p2c1_w(j),:);    
    end

    p1c2_layup_w = cell(length(mapped_p1c2_w),4);
    for j = 1:length(mapped_p1c2_w)
        p1c2_layup_w(j,:) = parent1.layup_w(mapped_p1c2_w(j),:);    
    end

    p2c2_layup_w = cell(length(mapped_p2c2_w),4);
    for j = 1:length(mapped_p2c2_w)
        p2c2_layup_w(j,:) = parent2.layup_w(mapped_p2c2_w(j),:);    
    end


    %% 4. Create empty layup cell arrays
    c1_layup_f1 = cell(length(mapped_p1c1_f1),4);
    c2_layup_f1 = cell(length(mapped_p1c2_f1),4);
    c1_layup_f2 = cell(length(mapped_p1c1_f2),4);
    c2_layup_f2 = cell(length(mapped_p1c2_f2),4);
    c1_layup_w  = cell(length(mapped_p1c1_w),4);
    c2_layup_w  = cell(length(mapped_p1c2_w),4);

    %% 5. Uniform crossover for material name and properties (thickness, density)
    %% nplies_f1
    % Get crossover values for this parameter
    if child1.nplies_f1 >= child2.nplies_f1
        randvals = unifrnd(0,1,child1.nplies_f1,1);
    else
        randvals = unifrnd(0,1,child2.nplies_f1,1);
    end

    % Child 1
    for j = 1:child1.nplies_f1
        randval = randvals(j);
        if randval <= crossThres
            c1_layup_f1(j,1) = p2c1_layup_f1(j,1);
            c1_layup_f1(j,3:4) = p2c1_layup_f1(j,3:4);
        else
            c1_layup_f1(j,1) = p1c1_layup_f1(j,1);
            c1_layup_f1(j,3:4) = p1c1_layup_f1(j,3:4);
        end
    end

    % Child 2
    for j = 1:child2.nplies_f1
        randval = randvals(j);
        if randval <= crossThres
            c2_layup_f1(j,1) = p1c2_layup_f1(j,1);
            c2_layup_f1(j,3:4) = p1c2_layup_f1(j,3:4);
        else
            c2_layup_f1(j,1) = p2c2_layup_f1(j,1);
            c2_layup_f1(j,3:4) = p2c2_layup_f1(j,3:4);
        end
    end

    %% nplies_f2
    % Get crossover values for this parameter
    if child1.nplies_f2 >= child2.nplies_f2
        randvals = unifrnd(0,1,child1.nplies_f2,1);
    else
        randvals = unifrnd(0,1,child2.nplies_f2,1);
    end

    % Child 1
    for j = 1:child1.nplies_f2
        randval = randvals(j);
        if randval <= crossThres
            c1_layup_f2(j,1) = p2c1_layup_f2(j,1);
            c1_layup_f2(j,3:4) = p2c1_layup_f2(j,3:4);
        else
            c1_layup_f2(j,1) = p1c1_layup_f2(j,1);
            c1_layup_f2(j,3:4) = p1c1_layup_f2(j,3:4);
        end
    end

    % Child 2
    for j = 1:child2.nplies_f2
        randval = randvals(j);
        if randval <= crossThres
            c2_layup_f2(j,1) = p1c2_layup_f2(j,1);
            c2_layup_f2(j,3:4) = p1c2_layup_f2(j,3:4);
        else
            c2_layup_f2(j,1) = p2c2_layup_f2(j,1);
            c2_layup_f2(j,3:4) = p2c2_layup_f2(j,3:4);
        end
    end

    %% nplies_w
    % Get crossover values for this parameter
    if child1.nplies_w >= child2.nplies_w
        randvals = unifrnd(0,1,child1.nplies_w,1);
    else
        randvals = unifrnd(0,1,child2.nplies_w,1);
    end

    % Child 1
    for j = 1:child1.nplies_w
        randval = randvals(j);
        if randval <= crossThres
            c1_layup_w(j,1) = p2c1_layup_w(j,1);
            c1_layup_w(j,3:4) = p2c1_layup_w(j,3:4);
        else
            c1_layup_w(j,1) = p1c1_layup_w(j,1);
            c1_layup_w(j,3:4) = p1c1_layup_w(j,3:4);
        end
    end

    % Child 2
    for j = 1:child2.nplies_w
        randval = randvals(j);
        if randval <= crossThres
            c2_layup_w(j,1) = p1c2_layup_w(j,1);
            c2_layup_w(j,3:4) = p1c2_layup_w(j,3:4);
        else
            c2_layup_w(j,1) = p2c2_layup_w(j,1);
            c2_layup_w(j,3:4) = p2c2_layup_w(j,3:4);
        end
    end

    %% 6. Uniform crossover for ply angles
    %% nplies_f1
    % Get crossover values for this parameter
    if child1.nplies_f1 >= child2.nplies_f1
        randvals = unifrnd(0,1,child1.nplies_f1,1);
    else
        randvals = unifrnd(0,1,child2.nplies_f1,1);
    end

    % Child 1
    for j = 1:child1.nplies_f1/2
        randval = randvals(j);
        if randval <= crossThres
            c1_layup_f1(j,2) = p2c1_layup_f1(j,2);
        else
            c1_layup_f1(j,2) = p1c1_layup_f1(j,2);
        end
    end
    c1_layup_f1((child1.nplies_f1/2)+1:end,:) = flip(c1_layup_f1(1:child1.nplies_f1/2,:));

    % Child 2
    for j = 1:child2.nplies_f1/2
        randval = randvals(j);
        if randval <= crossThres
            c2_layup_f1(j,2) = p1c2_layup_f1(j,2);
        else
            c2_layup_f1(j,2) = p2c2_layup_f1(j,2);
        end
    end
    c2_layup_f1((child2.nplies_f1/2)+1:end,:) = flip(c2_layup_f1(1:child2.nplies_f1/2,:));

    %% nplies_f2
    % Get crossover values for this parameter
    if child1.nplies_f2 >= child2.nplies_f2
        randvals = unifrnd(0,1,child1.nplies_f2,1);
    else
        randvals = unifrnd(0,1,child2.nplies_f2,1);
    end

    % Child 1
    for j = 1:child1.nplies_f2/2
        randval = randvals(j);
        if randval <= crossThres
            c1_layup_f2(j,2) = p2c1_layup_f2(j,2);
        else
            c1_layup_f2(j,2) = p1c1_layup_f2(j,2);
        end
    end
    c1_layup_f2((child1.nplies_f2/2)+1:end,:) = flip(c1_layup_f2(1:child1.nplies_f2/2,:));

    % Child 2
    for j = 1:child2.nplies_f2/2
        randval = randvals(j);
        if randval <= crossThres
            c2_layup_f2(j,2) = p1c2_layup_f2(j,2);
        else
            c2_layup_f2(j,2) = p2c2_layup_f2(j,2);
        end
    end
    c2_layup_f2((child2.nplies_f2/2)+1:end,:) = flip(c2_layup_f2(1:child2.nplies_f2/2,:));

    %% nplies_w
    % Get crossover values for this parameter
    if child1.nplies_w >= child2.nplies_w
        randvals = unifrnd(0,1,child1.nplies_w,1);
    else
        randvals = unifrnd(0,1,child2.nplies_w,1);
    end

    % Child 1
    for j = 1:child1.nplies_w/2
        randval = randvals(j);
        if randval <= crossThres
            c1_layup_w(j,2) = p2c1_layup_w(j,2);
        else
            c1_layup_w(j,2) = p1c1_layup_w(j,2);
        end
    end
    c1_layup_w((child1.nplies_w/2)+1:end,:) = flip(c1_layup_w(1:child1.nplies_w/2,:));

    % Child 2
    for j = 1:child2.nplies_w/2
        randval = randvals(j);
        if randval <= crossThres
            c2_layup_w(j,2) = p1c2_layup_w(j,2);
        else
            c2_layup_w(j,2) = p2c2_layup_w(j,2);
        end
    end
    c2_layup_w((child2.nplies_w/2)+1:end,:) = flip(c2_layup_w(1:child2.nplies_w/2,:));


    %% 7. Apply nplies_same rule
    % Child 1
    halfsame = child1.nplies_same/2;

    c1_layup_f2(1:halfsame,:) = c1_layup_f1(1:halfsame,:);
    c1_layup_w(1:halfsame,:) = c1_layup_f1(1:halfsame,:);
    
    if child1.nplies_f2 == 2
        c1_layup_f2(2,:) = c1_layup_f2(1,:);
    else
        c1_layup_f2 = balance_layup(c1_layup_f2);
    end
    
    if child1.nplies_w == 2
        c1_layup_w(2,:) = c1_layup_w(1,:);
    else
        c1_layup_w = balance_layup(c1_layup_w);
    end

    % Child 2
    halfsame = child2.nplies_same/2;

    c2_layup_f2(1:halfsame,:) = c2_layup_f1(1:halfsame,:);
    c2_layup_w(1:halfsame,:) = c2_layup_f1(1:halfsame,:);
    
    if child2.nplies_f2 == 2
        c2_layup_f2(2,:) = c2_layup_f2(1,:);
    else
        c2_layup_f2 = balance_layup(c2_layup_f2);
    end
    
    if child2.nplies_w == 2
        c2_layup_w(2,:) = c2_layup_w(1,:);
    else
        c2_layup_w = balance_layup(c2_layup_w);
    end
    
end


function [symmetric_layup] = balance_layup(layup)
% Make a given layup symmetric
    symmetric_layup = layup;
    
    nsym = size(layup,1)/2;
    n = size(layup,1);
    
    symmetric_layup(nsym+1:n,:) = flip(layup(1:nsym,:));

end

function [A_f1,A_f2,A_w,A,ybar,d_f1,d_f2,d_w,I_f1,I_f2,I_w,I] = derive_properties(child)
    % Derive remaining property values after independent parameters have been
    % set by the preceding code.

    % Calculate cross-sectional area
    A_f1 = child.b_f1*child.t_f1;
    A_f2 = child.b_f2*child.t_f2;
    A_w = child.h_w*child.t_w;
    A = A_f1 + A_f2 + A_w;

    % Calculate position of neutral axis
    yf1 = child.t_f1/2;
    yf2 = child.t_f2/2;
    yw = child.h_w/2;

    sum = yf2*A_f2;
    sum = sum + (yf1 + child.h_w + child.t_f2)*A_f1;
    sum = sum + (yw + child.t_f2)*A_w;
    ybar = sum/A;

    % Calculate area moment of inertia
    d_f1 = yf1 + child.h_w + child.t_f2 - ybar;
    d_f2 = yf2 - ybar;
    d_w = yw + child.t_f2 - ybar;

    %Ix = (bh^3)/12
    I_f1 = (child.b_f1*child.t_f1^3)/12;
    I_f2 = (child.b_f2*child.t_f2^3)/12;
    I_w = (child.t_w*child.h_w^3)/12;

    I = I_f1 + A_f1*d_f1^2 + I_f2 + A_f2*d_f2^2 + I_w + A_w*d_w^2;

end