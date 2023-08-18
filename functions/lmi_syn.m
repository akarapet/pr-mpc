function [E,Y1,Y2] = lmi_syn(A1,A2,B,Q,R)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Qh = sqrt(Q); Rh = sqrt(R);

setlmis([])
e = lmivar(2,[2 2]); 
y1 = lmivar(2,[1 2]);
y2 = lmivar(2,[1 2]);


% 1st LMI 
lmiterm([-1 1 1 e],1,1);
lmiterm([-1 1 2 e],1,1); 
lmiterm([-1 1 3 -y1],1,1); 
lmiterm([-1 1 4 e],1,A1'); 
lmiterm([-1 1 4 -y1],1,B'); 
lmiterm([-1 2 2 0],inv(Q));
lmiterm([-1 3 3 0],inv(R));
lmiterm([-1 4 4 e],1,1);

% 2nd LMI 
lmiterm([-2 1 1 e],1,1);
lmiterm([-2 1 2 e],1,1); 
lmiterm([-2 1 3 -y2],1,1); 
lmiterm([-2 1 4 e],1,A2'); 
lmiterm([-2 1 4 -y2],1,B'); 
lmiterm([-2 2 2 0],inv(Q));
lmiterm([-2 3 3 0],inv(R));
lmiterm([-2 4 4 e],1,1);

% 3rd LMI
lmiterm([-3 1 1 e],1,1);
% lmiterm([3 1 1 0],0.005);


lmis = getlmis;
% c = mat2dec(lmis,eye(2));
% options = [1e-2,0,0,0,0];
% [~,xopt] = mincx(lmis,c,options);
% E = dec2mat(lmis,xopt,e);
% Y1 = dec2mat(lmis,xopt,y1);
% Y2 = dec2mat(lmis,xopt,y2);
[~,xfeas] = feasp(lmis);
E = dec2mat(lmis,xfeas,e);
Y1 = dec2mat(lmis,xfeas,y1);
Y2 = dec2mat(lmis,xfeas,y2);
end

