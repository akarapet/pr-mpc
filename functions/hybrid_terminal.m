clc; close all; clear all;

%% Parameters

A1 = [0.35 -0.6062; 0.6062 0.35];
A2 = [0.35 0.6062; -0.6062 0.35];
B = [0; 1];

Q = eye(2); R = 0.4;

bx = 5*ones(4,1); bu = 1;
A_state = [eye(2); -eye(2)];

%% K1 K2 given in the example, terminal invariant set is given 

K11 = [-0.611 -0.3572]; K12 = [0.611 -0.3572];
P1 = Lyap_fun(A1,A2,B,K11,K12,Q,R);
% [A1_ini,b1_ini] = safe_set(A_state, bx, bu, K11, K12);
% [A1_reach,b1_reach] = one_step_reach(A_state, bx, A1+B*K11, A2+B*K12, A1_ini,b1_ini);
[A1max,b1max] = max_reach(A_state, bx, bu, A1, A2, B, K11, K12);

%% K1 K2 given by Riccati equation
[~,K21,~] = idare(A1,B,Q,R); K21 = -K21;
[~,K22,~] = idare(A2,B,Q,R); K22 = -K22;
P2 = Lyap_fun(A1,A2,B,K21,K22,Q,R);

[A2max,b2max] = max_reach(A_state, bx, bu, A1, A2, B, K21, K22);

%% Synthesize K1 and K2 using LMI

[E,Y1,Y2] = lmi_syn(A1,A2,B,Q,R);
P3 = inv(E);
K31 = Y1*P3; K32 = Y2*P3;

[A3max,b3max] = max_reach(A_state, bx, bu, A1, A2, B, K31, K32);


%% K1 = K2 = [0 0] Terminal invariant set is the state constraint
K41 = [0 0]; K42 = K41;
P4 = Lyap_fun(A1,A2,B,K41,K42,Q,R);

[A4max,b4max] = max_reach(A_state, bx, bu, A1, A2, B, K41, K42);


%%
X0= [-4 4.6; 5.3 5; -4.5 2.7; -5 -5]';
x0 = X0(:,3);
cost = hybrid_cost(x0,A1,A2,B,Q,R,K41,K42,A_state,bx,bu);
