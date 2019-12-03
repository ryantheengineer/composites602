classdef beamdesign
    properties
        nplies_f1       % Number of plies in the upper flange
        nplies_f2       % Number of plies in the lower flange
        nplies_w        % Number of plies in the web
        nplies_same     % Number of plies that are common between the flanges and the web
        layup_f1        % Cell array of materials and ply angles for the upper flange
        layup_f2        % Cell array of materials and ply angles for the lower flange
        layup_w         % Cell array of materials and ply angles for the web
        t_f1            % Thickness of the upper flange
        t_f2            % Thickness of the lower flange
        t_w             % Thickness of the web
        b_f1            % Width (base length) of the upper flange
        b_f2            % Width (base length) of the lower flange
        h_w             % Height of the web
        A_f1            % Area of upper flange
        A_f2            % Area of lower flange
        A_w             % Area of web
        A               % Area of the entire beam cross-section
        ybar            % Height from bottom of cross-section to neutral axis
        d_f1            % Distance from neutral axis to centroid of upper flange
        d_f2            % Distance from neutral axis to centroid of lower flange
        d_w             % Distance from neutral axis to centroid of web
        I_f1            % Area moment of inertia of upper flange
        I_f2            % Area moment of inertia of lower flange
        I_w             % Area moment of inertia of web
        I               % Beam cross-section area moment of inertia
        fitnesses       % Fitness score column vector for the design
    end
    
    methods
        function obj = beamdesign()
            load('material_properties.mat');
            
            % Set the number of plies to be used in the flanges and the web
            maxplies = 50;
            obj.nplies_f1 = random_nplies(maxplies);
            obj.nplies_f2 = random_nplies(maxplies);
            obj.nplies_w = random_nplies(maxplies);
            obj.nplies_same = random_nplies(min([obj.nplies_f1 obj.nplies_f2 obj.nplies_w]));
            
            % Create initial layups
            obj.layup_f1 = random_layup(obj.nplies_f1,'material_properties.mat');
            obj.layup_f2 = random_layup(obj.nplies_f2,'material_properties.mat');
            obj.layup_w = random_layup(obj.nplies_w,'material_properties.mat');
            
            % Force the outer nplies_same layers to be the same across all
            % members, using layup_f1 as the model layup
            halfsame = obj.nplies_same/2;
            
            obj.layup_f2(1:halfsame,:) = obj.layup_f1(1:halfsame,:);
            obj.layup_w(1:halfsame,:) = obj.layup_f1(1:halfsame,:);
            
            f1_start = length(obj.layup_f1) - halfsame;
            f2_start = length(obj.layup_f2) - halfsame;
            w_start = length(obj.layup_w) - halfsame;
            obj.layup_f2(f2_start:end,:) = obj.layup_f1(f1_start:end,:);
            obj.layup_w(w_start:end,:) = obj.layup_f1(f1_start:end,:);
            
            % Add up the total thickness of each member
            obj.t_f1 = 0;
            obj.t_f2 = 0;
            obj.t_w = 0;
            
            for i = 1:obj.nplies_f1
                rowname = MaterialProperties.Properties.RowNames(obj.layup_f1(i,1));
                obj.t_f1 = obj.t_f1 + MaterialProperties.t(rowname);                
            end
            
            for i = 1:obj.nplies_f2
                rowname = MaterialProperties.Properties.RowNames(obj.layup_f2(i,1));
                obj.t_f2 = obj.t_f2 + MaterialProperties.t(rowname);                
            end
            
            for i = 1:obj.nplies_w
                rowname = MaterialProperties.Properties.RowNames(obj.layup_w(i,1));
                obj.t_w = obj.t_w + MaterialProperties.t(rowname);                
            end
            
            % Generate continuous values for b_f1, b_f2, and h_w subject to
            % the maximum cross-section envelope and the previously
            % determined member thicknesses
            obj.b_f1 = unifrnd(obj.t_w+0.001,0.5);
            obj.b_f2 = unifrnd(obj.t_w+0.001,0.5);
            obj.h_w = unifrnd(0.001,(0.2667-obj.t_f1-obj.t_f2));
            
            % Calculate cross-sectional area
            obj.A_f1 = obj.b_f1*obj.t_f1;
            obj.A_f2 = obj.b_f2*obj.t_f2;
            obj.A_w = obj.h_w*obj.t_w;
            obj.A = obj.A_f1 + obj.A_f2 + obj.A_w;
            
            % Calculate position of neutral axis
            yf1 = obj.t_f1/2;
            yf2 = obj.t_f2/2;
            yw = obj.h_w/2;
            
            sum = yf2*obj.A_f2;
            sum = sum + (yf1 + obj.h_w + obj.t_f2)*obj.A_f1;
            sum = sum + (yw + obj.t_f2)*obj.A_w;
            obj.ybar = sum/obj.A;
            
            % Calculate area moment of inertia
            obj.d_f1 = yf1 + obj.h_w + obj.t_f2 - obj.ybar;
            obj.d_f2 = yf2 - obj.ybar;
            obj.d_w = yw + obj.t_f2 - obj.ybar;
            
            %Ix = (bh^3)/12
            obj.I_f1 = (obj.b_f1*obj.t_f1^3)/12;
            obj.I_f2 = (obj.b_f2*obj.t_f2^3)/12;
            obj.I_w = (obj.t_w*obj.h_w^3)/12;
            
            obj.I = obj.I_f1 + obj.A_f1*obj.d_f1^2 + obj.I_f2 + obj.A_f2*obj.d_f2^2 + obj.I_w + obj.A_w*obj.d_w^2;
            
            % Fitness of the multiple objectives
            obj.fitnesses = getFitness(obj);
        end
    end
end

%%% LOCAL FUNCTIONS %%%
function [nplies] = random_nplies(maxplies)
% Generate a random even integer number of plies (assumes we will
% always want a symmetric layup to follow design rules)
    while true
        % Generate a random integer
        nplies = randi([2 maxplies],1);
        
        % Check if the integer is even. If it is, exit the loop. Otherwise,
        % generate random integers until an even one is arrived at.
        if mod(nplies,2) == 0
            break
        else
            continue
        end
    end
end


function [layup] = random_layup(nplies,Properties_mat)
% Generate a random angle from an allowed list.
    
    % Load MaterialProperties table
    load(Properties_mat);
    
    layup = cell(nplies,4);
    
    % Define the allowable ply angles
    allowed_angles = [0 15 30 45 60 75 90 -15 -30 -45 -60 -75]';

    % Create a symmetric vector of ply angle layup
    angles = zeros(nplies,1);
    nsym = nplies/2;
    symplies = datasample(allowed_angles,nsym);
    for i = 1:nsym
        layup(i,2) = {symplies(i)};
    end
    
    symplies = flip(symplies);
    for i = (nsym+1):nplies
        layup(i,2) = {symplies(i-nsym)};
    end
    
    % Create a symmetric material layup
    symmaterials = datasample(MaterialProperties.Properties.RowNames,nsym);
    for i = 1:nsym
        layup(i,1) = symmaterials(i);
        rowname = MaterialProperties.Properties.RowNames(layup(i,1));
        layup(i,3) = {MaterialProperties.t(rowname)};
        layup(i,4) = {MaterialProperties.rho(rowname)};
    end
    
    symmaterials = flip(symmaterials);
    for i = (nsym+1):nplies
        layup(i,1) = symmaterials(i-nsym);
        rowname = MaterialProperties.Properties.RowNames(layup(i,1));
        layup(i,3) = {MaterialProperties.t(rowname)};
        layup(i,4) = {MaterialProperties.rho(rowname)};
    end
    
end


function [fitnesses] = getFitness(MaterialProperties,layup_f1,layup_f2,layup_w,t_f1,t_f2,t_w,b_f1,b_f2,h_w,delta_f1,delta_f2,delta_w,I_f1,I_f2,I_w,A_f1,A_f2,A_w,d_f1,d_f2,d_w,I)
% Fitness scoring function to determine how good a given design is, subject
% to the design constraints. Outputs fitnesses, a vector of the individual
% fitness measures of the objectives. This output is for a single design,
% which can be fed into a maximin fitness function later for a more
% comprehensive 

load('material_properties.mat');

%% Objective 1: Factor of safety above maximum bending load
FSreq = 1.5; % Required factor of safety

% Upper flange bending stiffness
[A_f1,B_f1,D_f1,alpha_f1,beta_f1,delta_f1] = compute_matrices(layup_f1,MaterialProperties);
EI_f1 = get_bend_stiff(t_f1,delta_f1,I_f1,A_f1,d_f1);

% Lower flange bending stiffness
[A_f2,B_f2,D_f2,alpha_f2,beta_f2,delta_f2] = compute_matrices(layup_f2,MaterialProperties);
EI_f2 = get_bend_stiff(t_f2,delta_f2,I_f2,A_f2,d_f2);

% Web bending stiffness
[A_w,B_w,D_w,alpha_w,beta_w,delta_w] = compute_matrices(layup_w,MaterialProperties);
EI_w = get_bend_stiff(t_w,delta_w,I_w,A_w,d_w);

% Equivalent bending stiffness
E_eq_bend = (EI_f1 + EI_f2 + EI_w)/beamdesign.I; % CHECK THIS
c_f1 = beamdesign.t_f2 + beamdesign.h_w + beamdesign.t_f1 - beamdesign.ybar;
c_f2 = beamdesign.ybar;
c = max(c_f1,c_f2);

maxmoment = max_bending_moment(E_eq_bend,c,beamdesign.I); % Calculate the maximum bending moment for the cross-section
moment_objective = -maxmoment; % Make moment_objective negative for maximin fitness? Not sure if this is necessary

%% Objective 2: Minimize weight
weight_objective = get_weight(beamdesign);

%% Objective 3: Minimize deflection (must be below 1 inch deflection)


%% Objective 4: Crippling


fitnesses = [moment_objective; weight_objective];
end

%%%%%%%%%% getFitness functions %%%%%%%%%%%
function [M] = max_bending_moment(sigma,c,I)
    M = sigma*I/c;
end


function [EImember] = get_bend_stiff(t,delta,Imember,Amember,dmember)
    Eb = 12/((t^3)*delta(1,1));
    EImember = Eb*(Imember + Amember*dmember^2);
end

function [section_weight] = get_weight(nplies_f1,nplies_f2,nplies_w,b_f1,b_f2,h_w,layup_f1,layup_f2,layup_w,A_f1,A_f2,A_w)
    % Compute the total mass per unit length of the beam
    
    % Upper flange
    w_f1 = 0;
    for i = 1:nplies_f1
        Aply = b_f1*cell2mat(layup_f1(i,3));
        w_f1 = w_f1 + Aply*cell2mat(layup_f1(i,4));
    end
    w_f1 = w_f1/A_f1;
    
    % Lower flange
    w_f2 = 0;
    for i = 1:nplies_f2
        Aply = b_f2*cell2mat(layup_f2(i,3));
        w_f2 = w_f2 + Aply*cell2mat(layup_f2(i,4));
    end
    w_f2 = w_f2/A_f2;
    
    % Web
    w_w = 0;
    for i = 1:nplies_w
        Aply = h_w*cell2mat(layup_w(i,3));
        w_w = w_w + Aply*cell2mat(layup_w(i,4));
    end
    w_w = w_w/A_w;
    
    % Total
    section_weight = w_f1 + w_f2 + w_w;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%