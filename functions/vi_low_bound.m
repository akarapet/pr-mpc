function vi_value = vi_low_bound(x0,A1,B1,A2,B2,Q,R,x_max,x_min,u_max,u_min,iteration)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
lti1 = LTISystem('A', A1, 'B', B1);
lti2 = LTISystem('A', A2, 'B', B2);

% choosing the discrete mode

A_u1 = [1 0];
A_u2 = [-1 0];

b_u = 0;

R_1 = Polyhedron('A',A_u1,'b',b_u);
R_2 = Polyhedron('A',A_u2,'b',b_u);


lti1.setDomain('u',R_1);
lti2.setDomain('u',R_2)

pwa = PWASystem([lti1,lti2]);

g = 2^iteration;
N = iteration;

for i = 1:g
    
    ctrls(i) = MPCController(pwa,N);
    % Add constraints on predicted states
    ctrls(i).model.x.min = [x_min;x_min];
    ctrls(i).model.x.max = [x_max;x_max];
    ctrls(i).model.u.min = [u_min;u_min];
    ctrls(i).model.u.max = [u_max;u_max];

    % Use quadratic state penalty with identity weighting matrix
    ctrls(i).model.x.penalty = QuadFunction(Q);

    % Set quadratic input penalty with identity weighting matrix
    ctrls(i).model.u.penalty = QuadFunction(R);    
    
    binary_i = int2bit(i-1,iteration);
    Y = ctrls(i).toYALMIP();
    for j = 1:iteration
        if binary_i(j) == 0
            Y.constraints = [Y.constraints, -1 <= Y.variables.u(1, j) <= -1 ];
        else
            Y.constraints = [Y.constraints, 1 <= Y.variables.u(1, j) <= 1 ];
        end
    end
    ctrls(i).fromYALMIP(Y);
    
    end

    %loop{i} = ClosedLoop(ctrls(i),pwa);
costs = [];
for i = 1:g
    [~, ~, openloop] = ctrls(i).evaluate(x0);
    cost = openloop.cost;
    costs = [costs,cost];
end
vi_value = min(costs);
vi_value = vi_value - iteration*R(1);
end


