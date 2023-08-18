close all;
clear;
addpath('functions')
%% Design Parameters

% Prediction horizon
ell = 2;

% Control horizon
T = 50;

% Input cost matrix
R = 0.4;

% State cost matrix
Q = eye(2);


% number of stable controllers

g = 4;

%% System Setup

% System Model
A1 = [0.35 -0.6062; 0.6062 0.35];
A2 = [0.35 0.6062; -0.6062 0.35];

B = [0; 1];

lti1 = LTISystem('A', A1, 'B', B);
lti2 = LTISystem('A', A2, 'B', B);

A_R_1 = [-1 0];
A_R_2 = [1 0];

b_R_1 = 0;
b_R_2 = 0;

R_1 = Polyhedron('A',A_R_1,'b',b_R_1);
R_2 = Polyhedron('A',A_R_2,'b',b_R_2);

lti1.setDomain('x',R_1);
lti2.setDomain('x',R_2)

pwa = PWASystem([lti1,lti2]);

% LQR solution
[~,K21,~] = idare(A1,B,Q,R); 
[~,K22,~] = idare(A2,B,Q,R);

% Lyapunov Solution

K11 = [0.611 0.3572]; K12 = [-0.611 0.3572];

% Synthesize the controller using LMI

[E,Y1,Y2] = lmi_syn(A1,A2,B,Q,R);
P3 = inv(E);
K31 = -Y1*P3; K32 = -Y2*P3;

%P1 = Lyap_fun(A1,A2,B,-K11,-K12,Q,R);

% Feedback gains

K_m{1}(1,:) = K21;
K_m{1}(2,:) = K22;

K_m{2}(1,:) = K11;
K_m{2}(2,:) = K12;

K_m{3}(1,:) = K31;
K_m{3}(2,:) = K32;

K_m{4}(1,:) = [0 0];
K_m{4}(2,:) = [0 0];

%% Input Constraints polyhedron 

A_2 = [1 ; -1];
b_2 = [1;1];

P_input_constraints = Polyhedron('A',A_2,'b',b_2);

%% State Constraints polyhedron
A_3 = [1 0; 0 1];
b_3 = [5;5];

P_state_constraints = Polyhedron('A',[A_3;-1*A_3],'b',[b_3;b_3]);


%% PWA Parameters

bx = 5*ones(4,1); bu = 1;
A_state = [eye(2); -eye(2)];


%% System evolution with suboptimal terminal set
for i = 1:g
    
    ctrls(i) = MPCController(pwa,ell);

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
    P_Ly{i} = Lyap_fun(A1,A2,B,-K_m{i}(1,:),-K_m{i}(2,:),Q,R);
    Lyapunov_Penalty{i} = QuadFunction(P_Ly{i});
    

    % Set constraint

    [A_Ly,b_Ly] = max_reach(A_state, bx, bu, A1, A2, B, -K_m{i}(1,:), -K_m{i}(2,:));

    InvSet_Lyapunov{i} = Polyhedron('A', A_Ly, 'b', b_Ly);

    % add a terminal set constraint (see help SystemSignal/filter_terminalSet)  
    ctrls(i).model.x.with('terminalSet');
    ctrls(i).model.x.terminalSet = InvSet_Lyapunov{i};
    
    % add terminal constraint
    ctrls(i).model.x.with('terminalPenalty');
    ctrls(i).model.x.terminalPenalty = Lyapunov_Penalty{i};

    loop{i} = ClosedLoop(ctrls(i),pwa);
end


%% Closed Loop Simulation For a grid of Initial condition


%Initial States
x1 = [-4, 5, -4.5, -5];
x2 = [4.6, 3.4, 2.7, -5];


X = [x1', x2'];

% closed loop cost of PRMPC
J_tilde_mu_M  = [];
J_mu_4 = [];

% MPC cost of controllers and T^lbar(J)
tilde_J_M = zeros(size(x1,2),g+1);



for i = 1:size(x1,2)

    x = [x1(i);x2(i)];
    
    [J_tilde_mu,tilde_J] = closed_loop_cost_calculator(x,T,g,ctrls,Q,R,loop);

    % closed loop cost
    J_tilde_mu_M = [J_tilde_mu_M;J_tilde_mu(g+1)];
   
    % MPC costs
    for j = 1:g+1
        tilde_J_M(i,j)= tilde_J(j);
    end

    % Calculate J_{mu_4}
    J_mu_4 = [J_mu_4; hybrid_cost(x,A1,A2,B,Q,R,K_m{g}(1,:),K_m{g}(2,:),A_state,bx,bu)];
end



tbl = [tilde_J_M, J_tilde_mu_M, J_mu_4];

tbl = tbl(:,[1 2 3 4 7 5 6]);
disp(round(tbl,1));
%% Functions

% J contains the closed loop cost in the following order [lqr,k1,...kn,multi-mpc]

function [J,tilde_J]= closed_loop_cost_calculator(x1,T,g,ctrls,Q,R,loop)



% MPC cost of controllers and the last one being T^l\bar{J}
tilde_J = ones(g+1,1)*Inf;

for j =1:g+1
    x1_m(:,j) = x1;
end

%  closed loop cost of PRMPC
J = ones(g+1,1)*0;


u = zeros(g,1);
% Closed Loop Evolution
for i = 1:T-1
    
    for j = 1:g

    
            [u(j),feas,~] = ctrls(j).evaluate(x1_m(:,j));
            
            [u_m(j),~,mpc_costs(j)] = ctrls(j).evaluate(x1_m(:,g+1));
            costs(j) = mpc_costs(j).cost;
            if i == 1
                tilde_J(j) = costs(j);
            end
            if feas
                data = loop{j}.simulate(x1_m(:,j),1);
                u(j) = data.U;
                J(j) = J(j)+ x1_m(:,j)'*Q*x1_m(:,j) + u(j)'*R*u(j);
                x1_m(:,j) = data.X(:,2);%A*x1_m(:,j)+ B*u(j);
            end
     
    end

    
    [~,I] = min(costs);
    if i ==1
        tilde_J(g+1) = min(costs);
    end
    u_mult = u_m(I); 
    J(g+1) = J(g+1)+ x1_m(:,g+1)'*Q*x1_m(:,g+1) + u_mult'*R*u_mult;
    data_m = loop{I}.simulate(x1_m(:,g+1),1);
    x1_m(:,g+1) = data_m.X(:,2);%A*x1_m(:,g+1)+ B*u_mult;
end

% cost at t = T
for j = 1:g+1
    J(j) = J(j)+ x1_m(:,j)'*Q*x1_m(:,j) ;
end

end