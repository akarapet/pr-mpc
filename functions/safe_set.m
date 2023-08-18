function [Aini,bini] = safe_set(A_state, bx, bu, K1, K2)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
Aini = [K1; -K1; K2; -K2; A_state];
bini = [bu*ones(4,1);bx];
end

