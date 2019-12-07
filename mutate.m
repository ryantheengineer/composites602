function [mutated_child] = mutate(child, currentGeneration, totalGen, beta, MaterialProperties)
% switchThres = 0.4; % Can change this
mutateThres = 0.2; % Can change this
alpha = getalpha(currentGeneration,totalGen,beta);

mutated_child = child;


%% 1. Mutate the number of nplies and round to the nearest even number
[child_nplies_f1,child_nplies_f2,child_nplies_w,child_nplies_same] = mutate_nplies(child,alpha);
mutated_child.nplies_f1 = child_nplies_f1;
mutated_child.nplies_f2 = child_nplies_f2;
mutated_child.nplies_w  = child_nplies_w;
mutated_child.nplies_same = child_nplies_same;


%% 2. Fix any discrepancies caused by changing nplies (like if nplies_same is too large or small)
min_nplies = min([mutated_child.nplies_f1 mutated_child.nplies_f2 mutated_child.nplies_w]);

if mutated_child.nplies_same > min_nplies
    mutated_child.nplies_same = min_nplies;
end


%% 3. If nplies increases, duplicate exterior plies of the layup. Otherwise, delete from the outside in.
diff_nplies_f1 = mutated_child.nplies_f1 - child.nplies_f1;
diff_nplies_f2 = mutated_child.nplies_f2 - child.nplies_f2;
diff_nplies_w = mutated_child.nplies_w - child.nplies_w;
diff_nplies_same = mutated_child.nplies_same - child.nplies_same;

temp_layup_f1 = cell(mutated_child.nplies_f1,4);
temp_layup_f2 = cell(mutated_child.nplies_f2,4);
temp_layup_w  = cell(mutated_child.nplies_w,4);

% layup_f1
if diff_nplies_f1 > 0
    % Place original child layup in the larger temp layup
    startind = diff_nplies_f1/2 + 1;
    endind = mutated_child.nplies_f1 - diff_nplies_f1/2;
    temp_layup_f1(startind:endind,:) = child.layup_f1;
    
    % Duplicate exterior plies
    for j = 1:diff_nplies_f1/2
        temp_layup_f1(j,:) = child.layup_f1(1,:);
    end
    for j = (mutated_child.nplies_f1 - diff_nplies_f1/2 + 1):mutated_child.nplies_f1
        temp_layup_f1(j,:) = child.layup_f1(1,:);
    end

elseif diff_nplies_f1 < 0
    startind = abs(diff_nplies_f1/2) + 1;
    endind = child.nplies_f1 + diff_nplies_f1/2;
    
    temp_layup_f1 = child.layup_f1(startind:endind,:);
else
    temp_layup_f1 = child.layup_f1;
end

mutated_child.layup_f1 = temp_layup_f1;

% layup_f2
if diff_nplies_f2 > 0
    % Place original child layup in the larger temp layup
    startind = diff_nplies_f2/2 + 1;
    endind = mutated_child.nplies_f2 - diff_nplies_f2/2;
    temp_layup_f2(startind:endind,:) = child.layup_f2;
    
    % Duplicate exterior plies
    for j = 1:diff_nplies_f2/2
        temp_layup_f2(j,:) = child.layup_f2(1,:);
    end
    for j = (mutated_child.nplies_f2 - diff_nplies_f2/2 + 1):mutated_child.nplies_f2
        temp_layup_f2(j,:) = child.layup_f2(1,:);
    end

elseif diff_nplies_f2 < 0
    startind = abs(diff_nplies_f2/2) + 1;
    endind = child.nplies_f2 + diff_nplies_f2/2;
    
    temp_layup_f2 = child.layup_f2(startind:endind,:);
else
    temp_layup_f2 = child.layup_f2;
end

mutated_child.layup_f2 = temp_layup_f2;

% layup_w
if diff_nplies_w > 0
    % Place original child layup in the larger temp layup
    startind = diff_nplies_w/2 + 1;
    endind = mutated_child.nplies_w - diff_nplies_w/2;
    temp_layup_w(startind:endind,:) = child.layup_w;
    
    % Duplicate exterior plies
    for j = 1:diff_nplies_w/2
        temp_layup_w(j,:) = child.layup_w(1,:);
    end
    for j = (mutated_child.nplies_w - diff_nplies_w/2 + 1):mutated_child.nplies_w
        temp_layup_w(j,:) = child.layup_w(1,:);
    end

elseif diff_nplies_w < 0
    startind = abs(diff_nplies_w/2) + 1;
    endind = child.nplies_w + diff_nplies_w/2;
    
    temp_layup_w = child.layup_w(startind:endind,:);
else
    temp_layup_w = child.layup_w;
end

mutated_child.layup_w = temp_layup_w;


%% 4. Mutate layup
[layup_f1,layup_f2,layup_w] = mutate_layup(mutated_child,MaterialProperties,alpha);
mutated_child.layup_f1 = layup_f1;
mutated_child.layup_f2 = layup_f2;
mutated_child.layup_w  = layup_w;


%% 5. Enforce nplies_same rule
halfsame = mutated_child.nplies_same/2;

mutated_child.layup_f2(1:halfsame,:) = mutated_child.layup_f1(1:halfsame,:);
mutated_child.layup_w(1:halfsame,:) = mutated_child.layup_f1(1:halfsame,:);


if mutated_child.nplies_f2 == 2
    mutated_child.layup_f2(2,:) = mutated_child.layup_f2(1,:);
else
    mutated_child.layup_f2 = balance_layup(mutated_child.layup_f2);
end

if mutated_child.nplies_w == 2
    mutated_child.layup_w(2,:) = mutated_child.layup_w(1,:);
else
    mutated_child.layup_w = balance_layup(mutated_child.layup_w);
end


%% 6. Mutate member length parameters
xmin = mutated_child.t_w + 0.001;
xmax = 0.5;
mutated_child.b_f1 = basicmutate(mutated_child.b_f1,xmin,xmax,alpha);
mutated_child.b_f2 = basicmutate(mutated_child.b_f2,xmin,xmax,alpha);

xmin = 0.001;
xmax = 0.2667 - mutated_child.t_f1 - mutated_child.t_f1;
mutated_child.h_w = basicmutate(mutated_child.h_w,xmin,xmax,alpha);


%% 7. Derive remaining parameter values
[t_f1,t_f2,t_w] = get_thicknesses(mutated_child,MaterialProperties);
mutated_child.t_f1 = t_f1;
mutated_child.t_f2 = t_f2;
mutated_child.t_w  = t_w;

[A_f1,A_f2,A_w,c_A,ybar,d_f1,d_f2,d_w,I_f1,I_f2,I_w,I] = derive_properties(mutated_child);
mutated_child.A_f1 = A_f1;
mutated_child.A_f2 = A_f2;
mutated_child.A_w  = A_w;
mutated_child.A    = c_A;
mutated_child.ybar = ybar;
mutated_child.d_f1 = d_f1;
mutated_child.d_f2 = d_f2;
mutated_child.d_w  = d_w;
mutated_child.I_f1 = I_f1;
mutated_child.I_f2 = I_f2;
mutated_child.I_w  = I_w;
mutated_child.I    = I;


end




%%%%%%%%%%%%%%%%%%%%%
%% Local functions %%
%%%%%%%%%%%%%%%%%%%%%
function [alpha] = getalpha(currentGeneration,totalGen,beta)
% Compute alpha value for dynamic mutation
    alpha = (1-(currentGeneration-1)/totalGen)^beta;
end

function [val] = basicmutate(x,xmin,xmax,alpha)
% Mutate the chosen parameter value
    r = unifrnd(xmin,xmax);    
    if r <= x
        val = xmin + ((r-xmin)^alpha)*((x-xmin)^(1-alpha));
    else
        val = xmax - ((xmax-r)^alpha)*((xmax-x)^(1-alpha));
    end
end

function [child_nplies_f1,child_nplies_f2,child_nplies_w,child_nplies_same] = mutate_nplies(child,alpha)
    % Mutate nplies parameter values    
    child_nplies_f1 = basicmutate(child.nplies_f1,2,50,alpha);
    child_nplies_f2 = basicmutate(child.nplies_f2,2,50,alpha);
    child_nplies_w  = basicmutate(child.nplies_w,2,50,alpha);
    
    nplies_set = [child_nplies_f1 child_nplies_f2 child_nplies_w];
    child_nplies_same = basicmutate(child.nplies_same,2,min(nplies_set),alpha);
    
    % Round nplies values to nearest even integer
    child_nplies_f1     = roundeven(child_nplies_f1);
    child_nplies_f2     = roundeven(child_nplies_f2);
    child_nplies_w      = roundeven(child_nplies_w);
    child_nplies_same   = roundeven(child_nplies_same);
    
end


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


function [layup_f1,layup_f2,layup_w] = mutate_layup(child,MaterialProperties,alpha)
% Put mutated_child into this function after the nplies parameters have
% been modified, so it only mutates plies that will remain in the final
% layup.
    layup_f1 = child.layup_f1;
    layup_f2 = child.layup_f2;
    layup_w  = child.layup_w;
    
    allowed_angles = [0 15 30 45 60 75 90 -15 -30 -45 -60 -75]';
    materials = MaterialProperties.Properties.RowNames;
    
    maxangle_ind = length(allowed_angles);
    minangle_ind = 1;
    
    maxmaterial_ind = size(MaterialProperties,1);
    minmaterial_ind = 1;
    
    %% Mutate materials
    % layup_f1
    for j = 1:child.nplies_f1/2
        xloc = find(ismember(materials,child.layup_f1(j,1)));
        newmatind = basicmutate(xloc,minmaterial_ind,maxmaterial_ind,alpha);
        newmatind = round(newmatind);
        layup_f1(j,1) = materials(newmatind);
        rowname = MaterialProperties.Properties.RowNames(layup_f1(j,1));
        layup_f1(j,3) = {MaterialProperties.t(rowname)};
        layup_f1(j,4) = {MaterialProperties.rho(rowname)};
    end
    
    % layup_f2
    for j = 1:child.nplies_f2/2
        xloc = find(ismember(materials,child.layup_f2(j,1)));
        newmatind = basicmutate(xloc,minmaterial_ind,maxmaterial_ind,alpha);
        newmatind = round(newmatind);
        layup_f2(j,1) = materials(newmatind);
        rowname = MaterialProperties.Properties.RowNames(layup_f2(j,1));
        layup_f2(j,3) = {MaterialProperties.t(rowname)};
        layup_f2(j,4) = {MaterialProperties.rho(rowname)};
    end
    
    % layup_w
    for j = 1:child.nplies_w/2
        xloc = find(ismember(materials,child.layup_w(j,1)));
        newmatind = basicmutate(xloc,minmaterial_ind,maxmaterial_ind,alpha);
        newmatind = round(newmatind);
        layup_w(j,1) = materials(newmatind);
        rowname = MaterialProperties.Properties.RowNames(layup_w(j,1));
        layup_w(j,3) = {MaterialProperties.t(rowname)};
        layup_w(j,4) = {MaterialProperties.rho(rowname)};
    end
    
    
    %% Mutate ply angles
    % layup_f1
    for j = 1:child.nplies_f1/2
        angle = cell2mat(child.layup_f1(j,2));
        xloc = find(allowed_angles == angle);
        newangleind = basicmutate(xloc,minangle_ind,maxangle_ind,alpha);
        newangleind = round(newangleind);
        layup_f1(j,2) = {allowed_angles(newangleind)};
    end
    
    % layup_f2
    for j = 1:child.nplies_f2/2
        angle = cell2mat(child.layup_f2(j,2));
        xloc = find(allowed_angles == angle);
        newangleind = basicmutate(xloc,minangle_ind,maxangle_ind,alpha);
        newangleind = round(newangleind);
        layup_f2(j,2) = {allowed_angles(newangleind)};
    end
    
    % layup_w
    for j = 1:child.nplies_w/2
        angle = cell2mat(child.layup_w(j,2));
        xloc = find(allowed_angles == angle);
        newangleind = basicmutate(xloc,minangle_ind,maxangle_ind,alpha);
        newangleind = round(newangleind);
        layup_w(j,2) = {allowed_angles(newangleind)};
    end
    
    %% Enforce symmetry
    if child.nplies_f1 == 2
        layup_f1(2,:) = layup_f1(1,:);
    else
        layup_f1 = balance_layup(layup_f1);
    end
    
    if child.nplies_f2 == 2
        layup_f2(2,:) = layup_f2(1,:);
    else
        layup_f2 = balance_layup(layup_f2);
    end
    
    if child.nplies_w == 2
        layup_w(2,:) = layup_w(2,:);
    else
        layup_w  = balance_layup(layup_w);
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

function [t_f1,t_f2,t_w] = get_thicknesses(child,MaterialProperties)
% Add up the total thickness of each member
    tply = MaterialProperties.t(1);
    
    t_f1 = child.nplies_f1*tply;
    t_f2 = child.nplies_f2*tply;
    t_w  = child.nplies_w*tply;

%     t_f1 = 0;
%     t_f2 = 0;
%     t_w = 0;
% 
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