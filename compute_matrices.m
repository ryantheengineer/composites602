function [A,B,D,tlaminate] = compute_matrices(layup,MaterialProperties)
% Compute the A, B, and D matrices for a given layup (2-dimensional matrix
% of material choices and fiber orientations) and material property table

nplies = size(layup,1);
tply = zeros(nplies,1);

R = [1 0 0; 0 1 0; 0 0 2];
A = zeros(3);
B = zeros(3);
D = zeros(3);

Qi = cell(nplies,1);

for i = 1:nplies
    % Select the appropriate material property
    material = layup(i,1);  % String corresponding to the RowName of the MaterialProperties table
    theta = layup(i,2);     % Ply angle in degrees
    
    % Fill out current basic material properties in x1-x2 coordinates
    E1 = MaterialProperties(material,'E1');
    E2 = MaterialProperties(material,'E2');
    G12 = MaterialProperties(material,'G12');
    nu12 = MaterialProperties(material,'nu12');
    nu23 = MaterialProperties(material,'nu23');
    t = MaterialProperties(material,'t');
    rho = MaterialProperties(material,'rho');
    
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
    
    T = [cosd(theta)^2,                      sind(theta)^2,       2*cosd(theta)*sind(theta);
         sind(theta)^2,                      cosd(theta)^2,      -2*cosd(theta)*sind(theta);
         -cosd(theta)*sind(theta), cosd(theta)*sind(theta), (cosd(theta)^2)-(sind(theta)^2)];
    
    Qi(i) = inv(T)*Qbasic*R*T*inv(R);
end

tlaminate = sum(tply);

z = zeros(nplies+1,1);

% Calculate the A, B, and D matrices
for i = 0:nplies
    z(i+1) = i*tply(i)-(tlaminate/2);
end

for k = 1:nplies-1
    A = A + Qi(k)*(z(k+1)-z(k));
    B = B + (1/2)*Qi(k)*(z(k+1)^2 - z(k)^2);
    D = D + (1/3)*Qi(k)*(z(k+1)^3 - z(k)^2);
end

end