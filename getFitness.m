function [fitnesses] = getFitness(design)
% Fitness scoring function to determine how good a given design is, subject
% to the design constraints. Outputs fitnesses, a vector of the individual
% fitness measures of the objectives. This output is for a single design,
% which can be fed into a maximin fitness function later for a more
% comprehensive 

load('material_properties.mat');

%% Objective 1: Factor of safety above maximum bending load
FSreq = 1.5; % Required factor of safety

% Upper flange bending stiffness
[~,~,~,~,~,delta_f1] = compute_matrices(design.layup_f1,MaterialProperties);
EI_f1 = get_bend_stiff(design.t_f1,delta_f1,design.I_f1,design.A_f1,design.d_f1);

% Lower flange bending stiffness
[~,~,~,~,~,delta_f2] = compute_matrices(design.layup_f2,MaterialProperties);
EI_f2 = get_bend_stiff(design.t_f2,delta_f2,design.I_f2,design.A_f2,design.d_f2);

% Web bending stiffness
[~,~,~,~,~,delta_w] = compute_matrices(design.layup_w,MaterialProperties);
EI_w = get_bend_stiff(design.t_w,delta_w,design.I_w,design.A_w,design.d_w);

% Equivalent bending stiffness
E_eq_bend = (EI_f1 + EI_f2 + EI_w)/design.I; % CHECK THIS
c_f1 = design.t_f2 + design.h_w + design.t_f1 - design.ybar;
c_f2 = design.ybar;
c = max(c_f1,c_f2);

maxmoment = max_bending_moment(E_eq_bend,c,design.I); % Calculate the maximum bending moment for the cross-section
allowmoment = 17460; % Max moment expected on the beam in normal use
FSmoment = allowmoment/maxmoment;
moment_objective = 1/FSmoment; % Make moment_objective negative for maximin fitness? Not sure if this is necessary

%% Objective 2: Minimize weight
weight_objective = get_weight(design.nplies_f1,design.nplies_f2,design.nplies_w,design.b_f1,design.b_f2,design.h_w,design.layup_f1,design.layup_f2,design.layup_w,design.A_f1,design.A_f2,design.A_w);

%% Objective 3: Minimize deflection (must be below 1 inch deflection)

%%LOAD CASE 8: (MOST EXTREME LOADING):
% load_8 =  [    97, -5383;
%             75.92, -3765;
%             40.55, -2626;
%             33.99, -1697;
%             11.33,  -251;
%                 0,   277;
%            -11.33,  1721;
%            -33.99,  2652;
%            -40.55,  3788;
%            -75.92,  5408;
%               -97, -5027];


%                 x pos     shear
load_8 = 1.0e+03*[      0,   5.027;
                  0.02108,   5.408;
                  0.05645,   3.788;
                  0.06301,   2.652;
                  0.08567,   1.721;
                    0.097,   0.277;
                  0.10833,  -0.251;
                  0.13099,  -1.697;
                  0.13755,  -2.626;
                  0.17292,  -3.765;
                    0.194,  -5.383];

length_scale_factor = .0254;
load_scale_factor = 4.44822;

p_vals = zeros(size(load_8,1)-1,2);
for i = 1:length(p_vals) %%this is where the center portion of the beam starts and ends
    p_vals(i,:) = [load_8(i+1,1)-load_8(1,1), load_8(i+1,2) - load_8(i,2)];
end

%%superimpose the p loads to find deflection
L = 2 * 97 * length_scale_factor; %total length from stanchion to stanchion

dmax = 0;
for i = 1:size(p_vals,1)
    b = L - (p_vals(i,1) * length_scale_factor); %distance from p to the right stanchion (gets smaller down the vector)
    single_dmax = (p_vals(i,2) * load_scale_factor)*b*(L^2 - b^2)^(3/2) / (9*sqrt(3) * E_eq_bend );
    dmax = dmax + single_dmax; %%superposition (conservative since there is also another downward force at the end which isn't being included)
end

m2in = 1/39.3701;
FS_deflection = m2in/abs(dmax);

dmax_objective = 1/FS_deflection;
% dmax_objective = FS_deflection;
% dmax_objective = dmax * 1.5;


%% Scale objectives to be on the same order of magnitude and combine into vector
moment_scale = 1000;
weight_scale = 5e2;
dmax_scale = 1e-4;

moment_objective = moment_objective/moment_scale;
weight_objective = weight_objective/weight_scale;
dmax_objective   = dmax_objective/dmax_scale;

fitnesses = [moment_objective; weight_objective; dmax_objective];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% getFitness functions %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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