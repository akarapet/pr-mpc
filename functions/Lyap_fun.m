function P = Lyap_fun(A1,A2,B,K1,K2,Q,R)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Ac1 = A1+B*K1; Ac2 = A2+B*K2; 
V1 = -Q-K1'*R*K1; V2 = -Q-K2'*R*K2; 
setlmis([])
p = lmivar(2,[2 2]); 

% 1st LMI 
lmiterm([-1 1 1 p],1,1);
lmiterm([-1 1 2 p],1,Ac1); 
lmiterm([-1 2 2 p],1,1);
lmiterm([-1 2 2 0],V1);

% 2nd LMI 
lmiterm([-2 1 1 p],1,1);
lmiterm([-2 1 2 p],1,Ac2);  
lmiterm([-2 2 2 p],1,1);
lmiterm([-2 2 2 0],V2);

% 3rd LMI
lmiterm([-3 1 1 p],1,1);

lmis = getlmis;
c = mat2dec(lmis,eye(2));
options = [1e-2,0,0,0,0];
[~,xopt] = mincx(lmis,c,options);
P = dec2mat(lmis,xopt,p);
end

