
% Solve for the voltages in an electric network for nonreversible Markov chains

syms R Q C l m lambda k real;

% Input voltage amplifiers as log-antisymmetric matrix:

lambda = [  [1      1/l     m       1]
            [l      1       1       1\m]
            [1/m    1       1       l]
            [1      m       1/l     1]];
        
% Input conductances as symmetric matrix:

C = [   [0      R       R*(l+1)*(m-1)/((l-1)*(m+1))       Q]
        [R      0       Q       R*(l+1)*(m-1)/((l-1)*(m+1))]
        [R*(l+1)*(m-1)/((l-1)*(m+1))      Q       0       R]
        [Q      R*(l+1)*(m-1)/((l-1)*(m+1))       R       0]];
    
n = length(C); 
u = sym('u',[1,n],'real'); 
U_div = sym('U',[1,n],'real');

% Set start and end voltages
u(1) = 1; u(n) = 0;

for k=2:n-1
    % set up expressions for voltage dividers (Proposition 3):
    U_div(k)=simplify(dot(u,C(k,:).*lambda(k,:)./(1+lambda(k,:)))/dot(C(k,:),1./(1+lambda(k,:))));
end

% solve voltage divider simultaneous equations
usol = struct2cell(solve( u(2:n-1) == U_div(2:n-1), u(2:n-1)));
for k = 1:n-2
    u(k+1)=usol{k};
end

% Currents
currents = sym('i',[n,n]);
for x = 1:n
    for y = 1:n
        currents(x,y) = simplify(2*C(x,y)/(1+lambda(y,x))*(lambda(y,x)*u(x)-u(y)));
    end
end

% Effective resistance
Reff = simplify((u(1)-u(n))/sum(currents(1,:)))