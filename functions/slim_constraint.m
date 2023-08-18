function [Aslim,bslim] = slim_constraint(Ain,bin)
%Removing redundent linear inequality constraints
%   Detailed explanation goes here
Atem = Ain; btem = bin; 
for i = size(Ain,1):-1:1
    f = -Ain(i,:); 
    A = Atem; A(i,:) = [];
    b = btem; b(i,:) = [];
    [~,fval,flag] = linprog(f,A,b);
    if flag > 0
    if -fval <= btem(i)
        Atem(i,:) = [];
        btem(i,:) = [];
    end
    end
end
Aslim = Atem; bslim = btem;

