function [A_reach,b_reach] = one_step_reach(A_state, bx, Ac1, Ac2, At,bt)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
A_reach = [At*Ac1; At*Ac2; At; A_state];
b_reach = [bt; bt; bt; bx];
[A_reach,b_reach] = slim_constraint(A_reach,b_reach);
end

