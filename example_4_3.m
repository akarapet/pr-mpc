close all;
clear;

%% Design Parameters

% Prediction horizon
ell = 3;

% Control horizon
T = 50;

% Input cost matrix
R = 1;

% State cost matrix
Q = eye(2);

% Plot the heatmap (plots if true)
plot_heatmap = false;


% number of stable controllers

g = 4;

%% System Setup

% System Model
A = [1 1; 0 1];
B = [1; 0.5];

lti = LTISystem('A', A, 'B', B);


% LQR solution
[K,P_ARE,~] = dlqr(A,B,Q,R);


% choose the ellipse
e = 3;

% Feedback gains

K_m = zeros(g,2);

K_m(1,:) = K;
K_m(2,:) = [0.1 1.2];
K_m(3,:) = [0.2 0.7];
K_m(4,:) = [0.3 0.8];

% old K 
%K_m(4,:) = [0.3 0.1];


% generate the last controller with an LMI

% [E,Y1] = lmi_syn(A,B,Q,R);
% P5 = inv(E);
% K_m(5,:) = -Y1*P5;

%% Input Constraints polyhedron 

A_2 = [1 ; -1];
b_2 = [1;1];

P_input_constraints = Polyhedron('A',A_2,'b',b_2);

%% State Constraints polyhedron
A_3 = [1 0; 0 1];
b_3 = [5;5];

P_state_constraints = Polyhedron('A',[A_3;-1*A_3],'b',[b_3;b_3]);


%% System evolution with suboptimal terminal set
for i = 1:g
    
    ctrls(i) = MPCController(lti,ell);

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


    % add terminal constraint
    ctrls(i).model.x.with('terminalPenalty');
    ctrls(i).model.x.terminalPenalty = Lyapunov_Penalty{i};


    % add a terminal set constraint (see help SystemSignal/filter_terminalSet)  
    ctrls(i).model.x.with('terminalSet');
    ctrls(i).model.x.terminalSet = InvSet_Lyapunov{i};


end


ctrls(e).model.x.terminalSet = ellipse{e};


%plot(InvSet_Lyapunov{4})
%hold on
%plot(ellipse{4})

InvSet_Lyapunov{e} =  ellipse{e};

%% Closed Loop Simulation For Single Initial condition

% Initial State
x1 = [-5; 2.7];


J = closed_loop_cost_calculator(x1,T,g,ctrls,Q,R,A,B);


%% Initial Feasible Region Calculation

for j = 1:g

    FSet_Lyapunov(j) = InvSet_Lyapunov{j};

    for i = 1:ell

        FSet_Lyapunov(j) = lti.reachableSet('X',FSet_Lyapunov(j),'U',P_input_constraints,'N',1,'direction', 'backward');
        FSet_Lyapunov(j) = FSet_Lyapunov(j).intersect(P_state_constraints).minHRep();

    end


end


U = Union(FSet_Lyapunov);  


%% Plot the terminal invariant sets
figure
colors = {[0.00000  0.78824  0.80000], [0.00000  0.78824  0.80000], [0.00000  0.78824  0.80000] ,[0.00000  0.78824  0.80000]};

% P_state_constraints.plot('Color','lightgray','linestyle','--','linewidth',1.5)
% hold on
U.plot('color',[0.00000  0.84314  1.00000],'linewidth',2,'edgecolor','b')
 
hold on

for i = 1:g

    InvSet_Lyapunov{g-i+1}.plot('Color',colors{i},'LineStyle',':','linewidth',2.5)
    hold on

end
grid on
xlabel('$x_1$','interpreter','latex','FontSize', 15,'FontWeight','bold')
ylabel('$x_2$','interpreter','latex','FontSize', 15,'FontWeight','bold')

legend('Support of $J_{\tilde{\mu}}$','','','','Support of $\bar{J}$','interpreter','latex','FontSize', 15)





%% Functions

% J contains the closed loop cost in the following order [lqr,k1,...kn,multi-mpc]

function J = closed_loop_cost_calculator(x1,T,g,ctrls,Q,R,A,B)


for j =1:g+1
    x1_m(:,j) = x1;
end

J = zeros(g+1,1);
u = zeros(g,1);
% Closed Loop Evolution
for i = 1:T-1

    for j = 1:g

        [u(j),~,~] = ctrls(j).evaluate(x1_m(:,j));
        [u_m(j),~,mpc_costs(j)] = ctrls(j).evaluate(x1_m(:,g+1));
        costs(j) = mpc_costs(j).cost;


        J(j) = J(j)+ x1_m(:,j)'*Q*x1_m(:,j) + u(j)'*R*u(j);
        x1_m(:,j) = A*x1_m(:,j)+ B*u(j);
    end
    
    [~,I] = min(costs);
    u_mult = u_m(I);  
    J(g+1) = J(g+1)+ x1_m(:,g+1)'*Q*x1_m(:,g+1) + u_mult'*R*u_mult;
    x1_m(:,g+1) = A*x1_m(:,g+1)+ B*u_mult;
end

% cost at t = T
for j = 1:g+1
    J(j) = J(j)+ x1_m(:,j)'*Q*x1_m(:,j) ;
end

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
