close all;
yalmip('clear')
clear;

%% Design Parameters

% Prediction horizon
N = 3;

% Control horizon
T = 50;

% Input cost matrix
R = 1;

% State cost matrix
Q = eye(2);

% number of stable controllers
g = 4;

%% System Setup

% System Model
A = [1 1; 0 1];
B = [1; 0.5];

% Number of states
nx = size(A,1); 

% Number of inputs
nu = size(B,2); 

lti = LTISystem('A', A, 'B', B);


% LQR solution
[K,P_ARE,~] = dlqr(A,B,Q,R);

% Feedback gains

K_m = zeros(g,2);

K_m(1,:) = K;
K_m(2,:) = [0.1 1.2];
K_m(3,:) = [0.2 0.7];
K_m(4,:) = [0.3 0.8];

% choose the index of an ellipsoidal constraint
e = 3;


%% Input Constraints polyhedron

A_2 = [1 ; -1];
b_2 = [1;1];

P_input_constraints = Polyhedron('A',A_2,'b',b_2);

%% State Constraints polyhedron
A_3 = [1 0; 0 1];
b_3 = [5;5];

P_state_constraints = Polyhedron('A',[A_3;-1*A_3],'b',[b_3;b_3]);

% union of pilyhedra 
P_nc = [];
%% System evolution with suboptimal terminal set
for i = 1:g

    ctrls(i) = MPCController(lti,N);

    % Add constraints on predicted states
    ctrls(i).model.x.min = [-5; -5];
    ctrls(i).model.x.max = [5; 5];

    % Add constraints on predicted control inputs
    ctrls(i).model.u.min = -1;
    ctrls(i).model.u.max = 1;

    % Use quadratic state penalty with identity weighting matrix
    ctrls(i).model.x.penalty = QuadFunction(Q);

    % Set quadratic input penalty with identity weighting matrix
    ctrls(i).model.u.penalty = QuadFunction(R);


    % Discrete Lyapunov Equation solution
    P_Ly{i} = dlyap((A-B*K_m(i,:))',Q+K_m(i,:)'*R*K_m(i,:));
    Lyapunov_Penalty{i} = QuadFunction(P_Ly{i});

    % Ellipse calculation

    alpha(i) = find_ellipsoid(P_Ly{i},-K_m(i,:),A_2,b_2,[A_3;-1*A_3],[b_3;b_3]);

    x_p{i} = sdpvar(2, 1);
    constraints_p{i} = x_p{i}'*P_Ly{i}*x_p{i} <= alpha(i);
    ellipse{i} = YSet(x_p{i}, constraints_p{i});
    pointlist = plot(constraints_p{i},x_p{i});
    pointlist = cell2mat(pointlist);
    pointlist = (unique(pointlist','rows','stable'))';
    pointlist = [pointlist pointlist(:,1)];

    ellipse{i} = Polyhedron(pointlist(:,1:2:end)');
    ellipse{i}.computeHRep;
    % Set constraint

    A_1 = [1 0; 0 1; -K_m(i,:)];
    A_Ly = [A_1; -1*A_1];
    b_1 = [5;5;1];
    b_Ly = [b_1;b_1];

    P_set_constraint{i} = Polyhedron('A', A_Ly, 'b', b_Ly);

    % Positively invariant set with defined with the stabilizing K_stable
    terminal_system{i} = LTISystem('A', A-B*K_m(i,:), 'B', [0; 0]);

    terminal_system{i}.x.with('setConstraint');
    terminal_system{i}.x.setConstraint = P_set_constraint{i};

    InvSet_Lyapunov{i} = terminal_system{i}.invariantSet();


    % add terminal penalty
    ctrls(i).model.x.with('terminalPenalty');
    ctrls(i).model.x.terminalPenalty = Lyapunov_Penalty{i};


    % add a terminal set constraint (see help SystemSignal/filter_terminalSet)
    ctrls(i).model.x.with('terminalSet');
    ctrls(i).model.x.terminalSet = InvSet_Lyapunov{i};
    
    % Terminal constraints polyhedra
    P_nc = [P_nc,InvSet_Lyapunov{i}];

end

% change e-th constraint set to an ellipse
ctrls(e).model.x.terminalSet = ellipse{e};
InvSet_Lyapunov{e} =  ellipse{e};


%% Binary Variable Naive implementation with YALMIP

% sdpvars for state and control inputs
u_naive = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x_naive = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

% g binary variables
for i = 1:g
    p(g) = binvar(1);
end

constraints = [];
objective = 0;

ops = sdpsettings('verbose',0);


% stage constraints and costs
for k = 1:N
    objective = objective + x_naive{k}'*Q*x_naive{k} + u_naive{k}'*R*u_naive{k};
    constraints = [constraints, x_naive{k+1} == A*x_naive{k}+B*u_naive{k}];
    constraints = [constraints, -1 <= u_naive{k}<= 1, -5<=x_naive{k+1}<=5];
end

objective_terminal = 0;

% terminal constraints and costs (with binary variables)
for i= 1:g
    objective_terminal = objective_terminal + p(i)*x_naive{N+1}'*P_Ly{i}*x_naive{N+1};
    constraints = [constraints, p(i)*(InvSet_Lyapunov{i}.A*x_naive{N+1}-InvSet_Lyapunov{i}.b) <= zeros(length(InvSet_Lyapunov{i}.A),1)];
end
constraints = [constraints, sum(p) == 1];

objective = objective + objective_terminal;

controller = optimizer(constraints,objective,ops,x_naive{1},{u_naive{1},p});

%% MPT Naive implementation with direct \bar{J}
ctrls_naive = MPCController(lti,N);
% Add constraints on predicted states
ctrls_naive.model.x.min = [-5; -5];
ctrls_naive.model.x.max = [5; 5];

% Add constraints on predicted control inputs
ctrls_naive.model.u.min = -1;
ctrls_naive.model.u.max = 1;

% Use quadratic state penalty with identity weighting matrix
ctrls_naive.model.x.penalty = QuadFunction(Q);

% Set quadratic input penalty with identity weighting matrix
ctrls_naive.model.u.penalty = QuadFunction(R);

% add a terminal set constraint (not necessary if J_bar is used)
%Y = ctrls_naive.toYALMIP();
%Y.constraints = Y.constraints + [ ismember( Y.variables.x(:, end), P_nc) ];
%Y.objective = Y.objective + term_func(Y.variables.x(:, end));
%ctrls_naive.fromYALMIP(Y);

% function J_bar
J_bar = @(x) min([max(1/double(~isnan(ismember(x,InvSet_Lyapunov{1})))-2,x'*P_Ly{1}*x),...
     max(1/double(~isnan(ismember(x,InvSet_Lyapunov{2})))-2,x'*P_Ly{2}*x),...
     max(1/double(~isnan(ismember(x,InvSet_Lyapunov{3})))-2,x'*P_Ly{3}*x),...
     max(1/double(~isnan(ismember(x,InvSet_Lyapunov{4})))-2,x'*P_Ly{4}*x)]);
F = Function(J_bar);

% add terminal penalty
ctrls_naive.model.x.with('terminalPenalty');
ctrls_naive.model.x.terminalPenalty = F;

%% Closed Loop Simulation

%Initial State
x1 = [-5.0 2.3; 2.7 -0.6];
sim_number = 100;

CL_costs=0;

average_trajectory_PR = [];
average_trajectory_naive = [];

for j = 1:length(x1)
    for i = 1 : sim_number
        [~,Time{i},Time_naive{i}] = closed_loop_cost_calculator(x1(:,j),T,g,ctrls,controller,Q,R,A,B);
        average_Time_PR(i,j) = mean(mean(Time{i}));
        average_Time_naive(i,j) = mean(Time_naive{i});

        average_trajectory_PR = [average_trajectory_PR,mean(Time{i},2)];
        average_trajectory_naive = [average_trajectory_naive, mean(Time_naive{i},2)];
    end
end

mean_PR = mean(average_Time_PR(2:end));
mean_naive = mean(average_Time_naive(2:end));

average_trajectory_PR_mean = mean(average_trajectory_PR(2:end,:),2);
average_trajectory_naive_mean = mean(average_trajectory_naive(2:end,:),2);

stdv1_curve_PR = average_trajectory_PR_mean + std(average_trajectory_PR(2:end,:),0,2);
stdv2_curve_PR = average_trajectory_PR_mean - std(average_trajectory_PR(2:end,:),0,2);

stdv1_curve_naive = average_trajectory_naive_mean + std(average_trajectory_naive(2:end,:),0,2);
stdv2_curve_naive = average_trajectory_naive_mean - std(average_trajectory_naive(2:end,:),0,2);

fprintf('Parallel Rollout MPC average computational time \n %f \n',mean(mean_PR));
fprintf('Naive MPC with binary variables average computational time \n %f \n',mean(mean_naive));

fprintf('PR-MPC is %f times faster than the naive implementation. \n',mean(mean_naive)/mean(mean_PR));
%% Plot
x_axis = 1:size(average_trajectory_PR_mean,1);

stdv_PR = [stdv1_curve_PR', flip(stdv2_curve_PR')];
fill ([x_axis, flip(x_axis)], stdv_PR, [0.3010 0.7450 0.9330],'FaceAlpha', 0.2)

hold on

plot(average_trajectory_PR_mean, LineWidth=3,Color="blue")

hold on

stdv_naive = [stdv1_curve_naive', flip(stdv2_curve_naive')];
fill ([x_axis, flip(x_axis)], stdv_naive, [0.6350 0.0780 0.1840],'FaceAlpha', 0.2)

hold on

plot(average_trajectory_naive_mean,LineWidth=3,Color="red")

legend("Standard deviation PR ", "Mean PR", "Standard deviation Naive", "Mean Naive", FontSize=20)
xlabel("$k$","Interpreter","latex",FontSize=20)
ylabel("Computational Time [seconds]", FontSize=20)


%% Functions

% J contains the closed loop cost in the following order [lqr,k1,...kn,pr-mpc, naive_multi_mpc]

function [J,Time, Time_naive] = closed_loop_cost_calculator(x1,T,g,ctrls,controller,Q,R,A,B)

% initial state
x1_naive = x1;

J_naive = 0;

u_naive = 0;

for j =1:g+1
    x1_m(:,j) = x1;
end

J = zeros(g+1,1);
u = zeros(g,1);
% Closed Loop Evolution

Time = zeros(T-1,1);
Time_naive = zeros(T-1,1);


for i = 1:T-1
    
    %tStart = tic;
    for j = 1:g
        tStart = tic;
        %[u(j),~,~] = ctrls(j).evaluate(x1_m(:,j));
        [u_m(j),~,mpc_costs(j)] = ctrls(j).evaluate(x1_m(:,g+1));
        costs(j) = mpc_costs(j).cost;
        %J(j) = J(j)+ x1_m(:,j)'*Q*x1_m(:,j) + u(j)'*R*u(j);
        %x1_m(:,j) = A*x1_m(:,j)+ B*u(j);
        Time(i,j)= toc(tStart);
    end
    %Time(i,:)= toc(tStart);
    [~,I] = min(costs);
    u_mult = u_m(I);
    
    

    J(g+1) = J(g+1)+ x1_m(:,g+1)'*Q*x1_m(:,g+1) + u_mult'*R*u_mult;
    x1_m(:,g+1) = A*x1_m(:,g+1)+ B*u_mult;


    % naive
    tStart_naive = tic;
    uk = controller{x1_naive};
    Time_naive(i) =  toc(tStart_naive);

    J_naive = J_naive + x1_naive'*Q*x1_naive + uk{1}'*R*uk{1};
    x1_naive = A*x1_naive + B*uk{1};

end

% cost at t = T
for j = 1:g+1
    J(j) = J(j)+ x1_m(:,j)'*Q*x1_m(:,j) ;
end

J_naive = J_naive+ x1_naive'*Q*x1_naive;

J = [J;J_naive];
end

function alpha = find_ellipsoid(P,K,H_u,h_u,H_x,h_x)


l_x = size(H_x,1);
l_u = size(H_u,1);

A =[];
b = [];

for i = 1:l_x

    A = [A;norm(inv(sqrtm(P))*H_x(i,:)')^2];
    b = [b;h_x(i)^2];
end

for i = 1:l_u

    A = [A;norm(inv(sqrtm(P))*K'*H_u(i,:)')^2];
    b = [b;h_u(i)^2];
end

alpha = linprog(-1,A,b);


end
