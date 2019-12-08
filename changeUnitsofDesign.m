function [t_f1,t_f2,t_w,b_f1,b_f2,h_w] = changeUnitsofDesign(design)
% Function for interpreting the design cross-section in inches instead of
% meters

m2in = 39.3701;

t_f1 = m2in*design.t_f1;
t_f2 = m2in*design.t_f2;
t_w  = m2in*design.t_w;

b_f1 = m2in*design.b_f1;
b_f2 = m2in*design.b_f2;
h_w  = m2in*design.h_w;


end