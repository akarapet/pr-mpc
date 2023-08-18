clc; clear all; close all;

%% Parameters
A = [1 1; 0 1]; B = [1; 0.5];
Q = eye(2); R = 1;

L1 = -[0.3 0.4]; L2 = -[1.5 0.2];

A1 = A+B*L1; A2 = A+B*L2;
Q1 = Q+L1'*R*L1; Q2 = Q+L2'*R*L2;
K1 = dlyap(A1',Q1); K2 = dlyap(A2',Q2);

K1_t = A'*K1*A-(A'*K1*B)*inv(R+B'*K1*B)*(B'*K1*A)+Q;
K2_t = A'*K2*A-(A'*K2*B)*inv(R+B'*K2*B)*(B'*K2*A)+Q;
DK = K1_t-K2_t;
[E,V]=eig(DK);
L1_t = -inv(R+B'*K1*B)*(B'*K1*A);
L2_t = -inv(R+B'*K2*B)*(B'*K2*A);

%% Plot
e = 0.001;
bound = 1.3;
[X,Y] = meshgrid(-bound:e:bound,-bound:e:bound);
Z = (DK(1)*X.^2+2*DK(2)*X.*Y+DK(4)*Y.^2>0)*1;
colormap summer
im=imagesc([X(1),X(end)],[Y(end),Y(1)],Z);
set(gca,'YDir','normal') 
im.AlphaData = .8;
xlabel('$x(1)$','interpreter','latex')
ylabel('$x(2)$','interpreter','latex')
hold on

N=20;
x = zeros(2,N); 
x(:,1)=[-0.8;0.9];
x1 = x; x2 = x;
for t=1:N
    x1(:,t+1) = A1*x1(:,t);
    x2(:,t+1) = A2*x2(:,t);
    
    if x(:,t)'*K1_t*x(:,t) < x(:,t)'*K2_t*x(:,t)
        x(:,t+1) = (A+ B*L1_t)*x(:,t);
    else
        x(:,t+1) = (A+ B*L2_t)*x(:,t);
    end
end
plot(x(1,:),x(2,:),'LineWidth',1.5);hold on;
plot(x1(1,:),x1(2,:),'--','LineWidth',1.5);hold on;
plot(x2(1,:),x2(2,:),':','LineWidth',1.5);hold on;
legend('Trajectory under $\tilde{\mu}$','Trajectory under $\mu_1$','Trajectory under $\mu_2$','interpreter','latex','AutoUpdate','off');

plot(x(1,1),x(2,1),'.','MarkerEdgeColor','blue','MarkerSize',15);hold on;
plot(x(1,2),x(2,2),'.','MarkerEdgeColor','blue','MarkerSize',15)