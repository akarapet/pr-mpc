close all;
clear;

%% Design Parameters

% Prediction horizon
N = 2;

% Control horizon
T = 12;

% Input cost matrix
R = 1;

% State cost matrix
Q = eye(2);

% Plot the heatmap (plots if true)
plot_heatmap = false;

% Heatmap gridding
delta = 0.01;

% number of stable controllers

g = 2;

%% System Setup

% System Model
A = [1 1; 0 1];
B = [1; 0.5];

lti = LTISystem('A', A, 'B', B);


% LQR solution
[K,P_ARE,e] = dlqr(A,B,Q,R);

% Feedback gains

K_m = zeros(g,2);

K_m(1,:) = K;
K_m(2,:) = [0.1 1.2];

% multiple prediction horizons

N = [2;3];


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
    
    ctrls(i) = MPCController(lti,N(i));

    % Add constraints on predicted states
    ctrls(i).model.x.min = [-5; -5];
    ctrls(i).model.x.max = [5; 5];

    % Add constraints on predicted control inputs
    ctrls(i).model.u.min = -1;
    ctrls(i).model.u.max = 1;

    % Use Inf-norm state penalty 
    ctrls(i).model.x.penalty = InfNormFunction(Q);

    % Set Inf-norm input penalty
    ctrls(i).model.u.penalty = InfNormFunction(R);    
    

    % Get the Final Cost to Go matrix
    %P_Ly{i} = dlyap((A-B*K_m(i,:))',Q+K_m(i,:)'*R*K_m(i,:));
    [~,P_Ly{i}] = P_matrix_1_inf(A-B*K_m(i,:),K_m(i,:),Q,R);
    Lyapunov_Penalty{i} = InfNormFunction(P_Ly{i});
    

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

    % add a terminal set constraint (see help SystemSignal/filter_terminalSet)  
    ctrls(i).model.x.with('terminalSet');
    ctrls(i).model.x.terminalSet = InvSet_Lyapunov{i};
    
    % add terminal constraint
    ctrls(i).model.x.with('terminalPenalty');
    ctrls(i).model.x.terminalPenalty = Lyapunov_Penalty{i};


end


%% Closed Loop Simulation For Single Initial condition


% Initial State
x1 = [-5; 2.7];


tilde_J_m = [];
Trajectory_cost = [];
UB = [];
for i = 1:T-1
    
    [J,x1_c,UB_e,I_mat,tilde_J] = closed_loop_cost_calculator(x1,T-i+1,g,ctrls,Q,R,A,B);
    Trajectory_cost =[Trajectory_cost;J(3)];
    UB = [UB;UB_e];
    tilde_J_m = [tilde_J_m;tilde_J'];
    x1= x1_c;
end
%%
figure 

plot(0:1:T-2,Trajectory_cost,'Color','#80B3FF','LineStyle','-','LineWidth',2)

hold on

plot(0:1:T-2,UB,'Color','#D95319','LineStyle','-','LineWidth',2)

xlabel('$k$','interpreter','latex','FontSize', 15)
ylabel('$J_{\tilde{\mu}}(x)$ $\quad$ / $\quad$ $T^\ell\bar{J}(x)$','interpreter','latex','FontSize', 15)

legend('$J_{\tilde{\mu}}(x)$','$T^\ell\bar{J}(x)$','interpreter','latex','FontSize', 15)

%% plot the bar graph

figure

bar(tilde_J_m)

xlabel('$k$','interpreter','latex','FontSize', 15)
ylabel('$T^{\ell_1}J_1(x)$ $\quad$ / $\quad$ $T^{\ell_2}J_2(x)$','interpreter','latex','FontSize', 15)

legend('$T^{\ell_1}J_1(x)$','$T^{\ell_2}J_2(x)$','interpreter','latex','FontSize', 15)

%% Initial Feasible Region Calculation

for j = 1:g

    FSet_Lyapunov(j) = InvSet_Lyapunov{j};

    for i = 1:N(j)

        FSet_Lyapunov(j) = lti.reachableSet('X',FSet_Lyapunov(j),'U',P_input_constraints,'N',1,'direction', 'backward');
        FSet_Lyapunov(j) = FSet_Lyapunov(j).intersect(P_state_constraints).minHRep();

    end


end

U = Union(FSet_Lyapunov);  

%% Cost Calculation over the feasible region

x1 = -5:delta:5;
x2 = x1;

L = 10/delta+1;

for i = 1:g+1
    Cost_matrix{i} = ones(L,L)*NaN;
end

%Cost_diff = ones(L,L)*-0.6;

if plot_heatmap
   for i = 1:L
        i
        for j = 1:L
            x = [x1(i);x2(j)];

            if U.contains(x)

                J = closed_loop_cost_calculator(x,T,g,ctrls,Q,R,A,B);

                for k = 1:g+1
                    Cost_matrix{k}(j,i) = J(k);
                end
            end

        end
    end


end


for i = 1:g

    Cost_diff{i} = Cost_matrix{i} -  Cost_matrix{g+1};
    Cost_diff{i} = 100 * Cost_diff{i}./Cost_matrix{g+1};
    Cost_diff{i}(isnan(Cost_diff{i})) = -0.6;
    Cost_diff{i}(5/delta+1,5/delta+1) = 0.01;
end

for i = 1:g+1
    Cost_matrix{i}(isnan(Cost_matrix{i})) = -0.6;
end


%% Upper Bound Calculations Over a Grid

for i = 1:g
    P_m{i} =  ones(L,L)*-4;
end

for i = 1:L
    for j = 1:L
        x = [x1(i);x2(j)];
        for k = 1:g
            if FSet_Lyapunov(k).contains(x)
                P_m{k}(j,i) = max(abs(P_Ly{k}*x));
            end
        end

    end
end




%% Functions

% J contains the closed loop cost in the following order [lqr,k1,...kn,multi-mpc]


function [J,x1_c,UB_e,I_mat,tilde_J] = closed_loop_cost_calculator(x1,T,g,ctrls,Q,R,A,B)

% MPC cost of controllers and the last one being T^l\bar{J}
tilde_J = ones(g,1)*Inf;

for j =1:g+1
    x1_m(:,j) = x1;
end

% which controller is used
I_mat = [];

J = zeros(g+1,1);
u = zeros(g,1);
% Closed Loop Evolution
for i = 1:T-1

    for j = 1:g

        [u(j),~,~] = ctrls(j).evaluate(x1_m(:,j));
        [u_m(j),~,mpc_costs(j)] = ctrls(j).evaluate(x1_m(:,g+1));
        costs(j) = mpc_costs(j).cost;
        if i == 1
            tilde_J(j) = costs(j);
        end

        J(j) = J(j)+  max(abs(Q*x1_m(:,j))) + max(abs(R*u(j)));
        x1_m(:,j) = A*x1_m(:,j)+ B*u(j);
    end
    
    [UB,I] = min(costs);

    u_mult = u_m(I);  
    I_mat = [I_mat;I];
    J(g+1) = J(g+1)+  max(abs(Q*x1_m(:,g+1))) +  max(abs(R*u_mult));
    x1_m(:,g+1) = A*x1_m(:,g+1)+ B*u_mult;

    if i == 1
        x1_c =  x1_m(:,g+1);
        UB_e = UB;
    end
end

% cost at t = T
for j = 1:g+1
    J(j) = J(j)+ max(abs(Q*x1_m(:,j))) ;
end


end


function [P_1,P_inf] = P_matrix_1_inf(F,K,Q,R)



[Av,Ad] = eig(F);
H1 = [real(Ad(1)), imag(Ad(1)); -imag(Ad(1)), real(Ad(1))];

[H1v,H1d] = eig(H1);
H = blkdiag(H1,zeros(2));
[Hv,Hd] = eig(H);

Gamma1 = [1+1i, 0; 0, 1-1i]; Gamma2 = zeros(2);
Gamma = [Gamma1; Gamma2];

Fac = [1, 1; 1i, -1i];
Fh = blkdiag(Fac,eye(2));

Pinf_tild = real(Hv*inv(Fh)*(Fh*Gamma*inv(Fac))*inv(Av*inv(Fac)));
Pinf_tild_left = inv(Pinf_tild.'*Pinf_tild)*Pinf_tild.';
sigma_inf = 1 - norm(H, "inf");
alpha_inf = norm(Q*Pinf_tild_left,"inf");
beta_inf = norm(R*K*Pinf_tild_left,"inf");

P_inf = (alpha_inf+beta_inf)/sigma_inf*Pinf_tild;

sigma_1 = 1 - norm(H, 1);
alpha_1 = norm(Q*Pinf_tild_left,1);
beta_1 = norm(R*K*Pinf_tild_left,1);

P_1 = (alpha_1+beta_1)/sigma_1*Pinf_tild;

end
