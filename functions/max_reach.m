function [Amax,bmax] = max_reach(A_state, bx, bu, A1, A2, B, K1, K2)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
[At,bt] = safe_set(A_state, bx, bu, K1, K2);
Ac1 = A1+B*K1; Ac2 = A2+B*K2;
flag = true;
while flag
    [A_reach,b_reach] = one_step_reach(A_state, bx, Ac1, Ac2, At,bt);
    if size(A_reach) == size(At)
    if norm(A_reach-At,Inf) < 0.001
       flag = false; 
    end
    end
    At = A_reach; bt = b_reach;
end
Amax = A_reach; bmax = b_reach;

