function [Ax,bx] = invariant_set(A, L, x_max, x_min, u_max, u_min)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
model = LTISystem('A', A);
n = length(A);
Ai = [eye(n);-eye(n); L; -L];
bi = [x_max*ones(n,1);-x_min*ones(n,1);u_max;-u_min];
P = Polyhedron('A', Ai, 'b', bi);
model.x.with('setConstraint');
model.x.setConstraint = P;
InvSet = model.invariantSet();
Ax = InvSet.A;
bx = InvSet.b;
end

