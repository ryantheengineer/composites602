function [fitness] = getFitness(beamdesign)
% Fitness scoring function (maximin to deal with multiple objectives) to
% determine how good a given design is, subject to the design constraints.

load('material_properties.mat');

%% Objective 1: Factor of safety above maximum bending load
FSreq = 1.5; % Required factor of safety

% Upper flange bending stiffness
[A_f1,B_f1,D_f1,alpha_f1,beta_f1,delta_f1] = compute_matrices(beamdesign.layup_f1,MaterialProperties);
EI_f1 = get_bend_stiff(beamdesign.t_f1,delta_f1,beamdesign.I_f1,beamdesign.A_f1,beamdesign.d_f1);

% Lower flange bending stiffness
[A_f2,B_f2,D_f2,alpha_f2,beta_f2,delta_f2] = compute_matrices(beamdesign.layup_f2,MaterialProperties);
EI_f2 = get_bend_stiff(beamdesign.t_f2,delta_f2,beamdesign.I_f2,beamdesign.A_f2,beamdesign.d_f2);

% Web bending stiffness
[A_w,B_w,D_w,alpha_w,beta_w,delta_w] = compute_matrices(beamdesign.layup_w,MaterialProperties);
EI_w = get_bend_stiff(beamdesign.t_w,delta_w,beamdesign.I_w,beamdesign.A_w,beamdesign.d_w);

% Equivalent bending stiffness
EI_eq_bend = (EI_f1 + EI_f2 + EI_w)/beamdesign.I;

% maxmoment = ; % Specified maximum moment ()

%% Objective 2: Minimize weight
section_weight = get_weight(beamdesign);
fitness = section_weight;

%% Objective 3: Minimize deflection (must be below 1 inch deflection)


%% Objective 4: Crippling



%% Maximin fitness function to combine the different objectives and create
% a Pareto front



end


%% LOCAL FUNCTIONS %%
function [moment] = max_bending_moment(sigma,c,I)
    M = sigma*I/c;
end

function [EImember] = get_bend_stiff(t,delta,Imember,Amember,dmember)
    Eb = 12/((t^3)*delta(1,1));
    EImember = Eb*(Imember + Amember*dmember^2);
end

function [section_weight] = get_weight(beamdesign)
    % Compute the total mass per unit length of the beam
    
    % Upper flange
    w_f1 = 0;
    for i = 1:beamdesign.nplies_f1
        Aply = beamdesign.b_f1*cell2mat(beamdesign.layup_f1(i,3));
        w_f1 = w_f1 + Aply*cell2mat(beamdesign.layup_f1(i,4));
    end
    w_f1 = w_f1/beamdesign.A_f1;
    
    % Lower flange
    w_f2 = 0;
    for i = 1:beamdesign.nplies_f2
        Aply = beamdesign.b_f2*cell2mat(beamdesign.layup_f2(i,3));
        w_f2 = w_f2 + Aply*cell2mat(beamdesign.layup_f2(i,4));
    end
    w_f2 = w_f2/beamdesign.A_f2;
    
    % Web
    w_w = 0;
    for i = 1:beamdesign.nplies_w
        Aply = beamdesign.h_w*cell2mat(beamdesign.layup_w(i,3));
        w_w = w_w + Aply*cell2mat(beamdesign.layup_w(i,4));
    end
    w_w = w_w/beamdesign.A_w;
    
    % Total
    section_weight = w_f1 + w_f2 + w_w;
    
end