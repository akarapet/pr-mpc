clc; close all;
clear;

%% Design Parameters

% Prediction horizon
ell = 5;

% Control horizon
T = 80;

% Input cost matrix
R = [0.00001 0; 0 1];
R_s = 1;

% State cost matrix
Q = eye(2);




%% System Setup

% System Model
A1 = [2 1; 0 1];
A2 = [2 1; 0 0.5];
n = length(A1);
B1 = [0 1; 0 1];
B2 = [0 1; 0 2];
m = size(B1,2);
B1_s = [1;1];
B2_s = [1;2];

% x_max = Inf; x_min = -Inf;
% u_max = Inf; u_min = -Inf;

x_max = 5; x_min = -5;
u_max = 1; u_min = -1;

lti1 = LTISystem('A', A1, 'B', B1);
lti2 = LTISystem('A', A2, 'B', B2);


% choosing the discrete mode
% if u_1<0 => lti 1
% if u_1>0 => lti 2

A_u1 = [1 0];
A_u2 = [-1 0];

b_u = 0;

R_1 = Polyhedron('A',A_u1,'b',b_u);
R_2 = Polyhedron('A',A_u2,'b',b_u);


lti1.setDomain('u',R_1);
lti2.setDomain('u',R_2)

pwa = PWASystem([lti1,lti2]);

% LQR solution
[~,K21,~] = idare(A1,B1_s,Q,R_s); 
[~,K22,~] = idare(A2,B2_s,Q,R_s);

% Lyapunov solution

zeta = 1;
[~,K11,~] = idare(A1,B1_s,Q,zeta*R_s); 
[~,K12,~] = idare(A2,B2_s,Q,zeta*R_s);


% Feedback gains

K_m{1}(1,:) = K21;
K_m{1}(2,:) = K22;

K_m{2}(1,:) = K11;
K_m{2}(2,:) = K12;

% Discrete Lyapunov Equation solution
P_Ly{1} = dlyap((A1-B1_s*K_m{2}(1,:))',Q+K_m{2}(1,:)'*R_s*K_m{2}(1,:));
P_Ly{2} = P_Ly{1};
P_Ly{3} = dlyap((A2-B2_s*K_m{2}(2,:))',Q+K_m{2}(2,:)'*R_s*K_m{2}(2,:));
P_Ly{4} = P_Ly{3};

[A_i{1},b_i{1}] = invariant_set((A1-B1_s*K_m{2}(1,:)), -K_m{2}(1,:), x_max, x_min, u_max, u_min);
A_i{2}= A_i{1}; b_i{2} = b_i{1};
[A_i{3},b_i{3}] = invariant_set((A2-B2_s*K_m{2}(2,:)), -K_m{2}(2,:), x_max, x_min, u_max, u_min);
A_i{4}= A_i{3}; b_i{4} = b_i{3};

% number of modes
g = 4;

%% System evolution with suboptimal terminal set
for i = 1:g
    
    ctrls(i) = MPCController(pwa,ell);
    % Add constraints on predicted states
    ctrls(i).model.x.min = [x_min;x_min];
    ctrls(i).model.x.max = [x_max;x_max];

    % Use quadratic state penalty with identity weighting matrix
    ctrls(i).model.x.penalty = QuadFunction(Q);

    % Set quadratic input penalty with identity weighting matrix
    ctrls(i).model.u.penalty = QuadFunction(R);    


    % add terminal constraint
    ctrls(i).model.x.with('terminalPenalty');
    Lyapunov_Penalty{i} = QuadFunction(P_Ly{i});
    ctrls(i).model.x.terminalPenalty = Lyapunov_Penalty{i};

    if i == 1 % mode 1: 1 1 ... 1 Terminal_Set_1
        Y1 = ctrls(i).toYALMIP();
        Y1.constraints = [Y1.constraints, -1 <= Y1.variables.u(1, 1:end) <= -1 ]; % force lti 1
        Y1.constraints = [Y1.constraints, u_min <= Y1.variables.u(2, 1:end) <= u_max ];
        Y1.constraints = [Y1.constraints, A_i{i}*Y1.variables.x(:,end) <= b_i{i} ];
        ctrls(i).fromYALMIP(Y1);
    elseif i == 2 % mode 2: 2 1 ... 1 Terminal_Set_1
        Y2 = ctrls(i).toYALMIP();
        Y2.constraints = [Y2.constraints, 1 <= Y2.variables.u(1, 1) <= 1 ]; % force lti 2
        Y2.constraints = [Y2.constraints, -1 <= Y2.variables.u(1, 2:end) <= -1 ];
        Y2.constraints = [Y2.constraints, u_min <= Y2.variables.u(2, 1:end) <= u_max ];
        Y2.constraints = [Y2.constraints, A_i{i}*Y2.variables.x(:,end) <= b_i{i} ];
        ctrls(i).fromYALMIP(Y2);
    elseif i ==3 % mode 3: 2 2 ... 2 Terminal_Set_2
        Y3 = ctrls(i).toYALMIP();
        Y3.constraints = [Y3.constraints, 1 <= Y3.variables.u(1, 1:end) <= 1 ]; % force lti 2
        Y3.constraints = [Y3.constraints, u_min <= Y3.variables.u(2, 1:end) <= u_max ];
        Y3.constraints = [Y3.constraints, A_i{i}*Y3.variables.x(:,end) <= b_i{i} ];
        ctrls(i).fromYALMIP(Y3);
    else % mode 4: 1 2 ... 2 Terminal_Set_2
        Y4 = ctrls(i).toYALMIP();
        Y4.constraints = Y4.constraints + [ -1 <= Y4.variables.u(1, 1) <= -1 ];% force lti 1
        Y4.constraints = Y4.constraints + [ 1 <= Y4.variables.u(1, 2:end) <= 1];
        Y4.constraints = Y4.constraints + [ u_min <= Y4.variables.u(2, 1:end) <= u_max];
        Y4.constraints = [Y4.constraints, A_i{i}*Y4.variables.x(:,end) <= b_i{i} ];
        ctrls(i).fromYALMIP(Y4);
    end

    % sloop{i} = ClosedLoop(ctrls(i),pwa);
end


%% openloop comparison

% initial state
x0 = [1;1];


% H{1} = [6.064 1.205; 1.205 1.905];
% H{2} = [9.084 3.233; 3.233 2.347];
% H{3} = [5.107 1.266; 1.266 1.935];
% H{4} = [7.216 2.560; 2.560 2.106];
% 
% Lower_bound = [];
% costs = [];
% for i = 1:g
% 
%     [~, ~, openloop] = ctrls(i).evaluate(x0);
%     cost = openloop.cost;
%     costs = [costs,cost];
%     Lower_bound(i) = x0'*H{i}*x0;
% 
% end
% 
% pr_cost_ol = min(costs);
% lb = min(Lower_bound);


%% closed loop

% number of VI iterations
iteration = 8;

%Initial States
x1 = [-4, 1.2, -3.5, -1.5];
x2 = [4.6, 1.5, 2.0, -0.5];
% x1 = [-4];
% x2 = [4.6]; 


X_0 = [x1', x2'];

% table for the switched system
tbl = [];

for k = 1:length(x1)

    X = X_0(k,:)'; U = [];
    J = 0;
    for j = 0:T
        x_current = X(:,end);
        X_tem = zeros(n,g); U_tem = zeros(m,g);
        costs_tem = zeros(1,g);
        for i = 1:g
            [~, ~, openloop] = ctrls(i).evaluate(x_current);
            X_tem(:,i) = openloop.X(:,2);
            U_tem(:,i) = openloop.U(:,1);
            costs_tem(i) = openloop.cost;
        end
        [~,index] = min(costs_tem);
        J = J + X(:,end)'*Q*X(:,end)+U_tem(:,index)'*R*U_tem(:,index);
        
        if j == 0
            tbl(k,1) = min(costs_tem);
        end

        X = [X,X_tem(:,index)];
        U = [U,U_tem(:,index)];
        J = J - R(1);
    end

    % Cost of PR_Algorithm
    J = J + X(:,end)'*P_Ly{1}*X(:,end);

    % Lower bound with constrained Value Iteration

    vi_value = vi_low_bound(X_0(k,:)',A1,B1,A2,B2,Q,R,x_max,x_min,u_max,u_min,iteration);

    % suboptimality gap

    sub_gap = (J - vi_value)/J;

    % Populate the table
    tbl(k,2) = J;
    tbl(k,3) = sub_gap;
    

end


fprintf('%3.1f %3.1f %1.1s \n',tbl')
