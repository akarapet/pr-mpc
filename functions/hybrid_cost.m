function cost = hybrid_cost(x0,A1,A2,B,Q,R,L1,L2,Ax,bx,bu)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
N=50;
x = zeros(2,N); 
x(:,1)= x0;
cost = 0;
for t=1:N
    if Ax*x(:,t) > bx
        cost = Inf; break;
    else
        if x(1,t)>= 0 
           if abs(L1*x(:,t))> bu
               cost = Inf; break;
           else
               cost = cost + x(:,t)'*(Q+L1'*R*L1)*x(:,t);
               x(:,t+1) = (A1+ B*L1)*x(:,t);
           end
        else
            if abs(L2*x(:,t))> bu
                cost = Inf; break;
            else
               cost = cost + x(:,t)'*(Q+L2'*R*L2)*x(:,t);
               x(:,t+1) = (A2+ B*L2)*x(:,t);
            end
        end
    end
end

end

