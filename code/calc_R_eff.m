function Reff = calc_R_eff(C, lambda)
% Calculate the effective resistance across a network. C = conductance
% matrix, lambda = amplifier matrix.

n = length(C);

% Set up vectors to contain voltage solutions and equations
u_sol = sym('u',[1,n],'real'); % voltage solutions
U_div = sym('U',[1,n],'real'); % voltage divider equations

% Set start and end voltages
u_sol(1) = 1; u_sol(n) = 0;

for k=2:n-1
    % set up expressions for voltage dividers:
    U_div(k)=simplify(dot(u_sol,C(k,:).*lambda(k,:)./(1+lambda(k,:)))/dot(C(k,:),1./(1+lambda(k,:))));
end

% solve voltage divider simultaneous equations to obtain voltages
usol = struct2cell(solve( u_sol(2:n-1) == U_div(2:n-1), u_sol(2:n-1)));
for k = 1:n-2
    u_sol(k+1)=usol{k};
end

% Currents
currents = sym('i',[n,n]);
for x = 1:n
    for y = 1:n
        currents(x,y) = simplify(2*C(x,y)/(1+lambda(y,x))*(lambda(y,x)*u_sol(x)-u_sol(y)));
    end
end

% Effective resistance
Reff = simplify((u_sol(1)-u_sol(n))/sum(currents(1,:)));