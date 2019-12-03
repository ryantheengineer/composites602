function [A,B,D,tlaminate] = compute_matrices(layup,MaterialProperties)
% Compute the A, B, and D matrices for a given layup (2-dimensional cell
% array of material choices and fiber orientations, with the first row
% corresponding to the bottom of the laminate) and material property table

nplies = size(layup,1);
tply = zeros(nplies,1);
material = layup(:,1); % String corresponding to the RowName of the MaterialProperties table
theta = cell2mat(layup(:,2)); % Ply angle in degrees

R = [1 0 0; 0 1 0; 0 0 2];
A = zeros(3);
B = zeros(3);
D = zeros(3);

Qi = zeros(3,3,nplies);

for i = 1:nplies
    % Fill out current basic material properties in x1-x2 coordinates
    E1 = MaterialProperties.E1(material(i));
    E2 = MaterialProperties.E2(material(i));
    G12 = MaterialProperties.G12(material(i));
    nu12 = MaterialProperties.nu12(material(i));
    nu23 = MaterialProperties.nu23(material(i));
    t = MaterialProperties.t(material(i));
    rho = MaterialProperties.rho(material(i));
    
    % Calculate the remaining material properties
    nu13 = nu12;
    nu21 = nu12;
    nu32 = nu23;
    G13 = G12;
    tply(i) = t;
    
    S = zeros(6);
    S(1:3,1:3) = [1/E1,     -nu21/E1, -nu21/E1;
                  -nu12/E1,     1/E2, -nu32/E2;
                  -nu12/E1, -nu23/E2,     1/E2];
    
    S(4:6,4:6) = [2*(1+nu23)/E2,    0,     0;
                  0,            1/G13,     0;
                  0,                0, 1/G13];
    
    
    Qbasic = inv([S(1,1), S(1,2), S(1,6);
                  S(1,2), S(2,2), S(2,6);
                  S(1,6), S(2,6), S(6,6)]);
    
    T = [cosd(theta(i))^2,                            sind(theta(i))^2,       2*cosd(theta(i))*sind(theta(i));
         sind(theta(i))^2,                            cosd(theta(i))^2,      -2*cosd(theta(i))*sind(theta(i));
         -cosd(theta(i))*sind(theta(i)), cosd(theta(i))*sind(theta(i)), (cosd(theta(i))^2)-(sind(theta(i))^2)];
    
    Qi(:,:,i) = inv(T)*Qbasic*R*T*inv(R);
end

tlaminate = sum(tply);

z = zeros(nplies+1,1);

% Calculate the A, B, and D matrices
for i = 2:length(z)
    z(i) = sum(tply(1:i-1))-sum(tply(1:(nplies/2))); % WORKS ONLY FOR SYMMETRIC LAMINATES
end

z(1) = z(1) - sum(tply(1:(nplies/2)));

for k = 1:nplies
    A = A + Qi(:,:,k)*(z(k+1)-z(k));
    B = B + (1/2)*Qi(:,:,k)*(z(k+1)^2 - z(k)^2);
    D = D + (1/3)*Qi(:,:,k)*(z(k+1)^3 - z(k)^3);
end

end