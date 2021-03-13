n=5;
trials=10;
R_eff = zeros(trials,2);
trial = 1;
while trial<=    trials %<- Number of networks to test

    %STEP 1
    %Generate matrix G_{xy} = D_{xy}*gamma_{xy} whose rows and columns sum to Const

    Const=1;

    G=rand(n).*not(eye(n));
    G=Const*bistochastic(G,0.01,200);

    % STEP 2
    % Calculate the gamma_{xy} and the initial D_{xy}

    D_init = zeros(n);
    gamma = ones(n);
    for i = 1:n
        for j = 1:n
            D_init(i,j) = sqrt(G(i,j)*G(j,i));
            if i~=j
                gamma(i,j) = sqrt(G(i,j)/G(j,i));
            end
        end
    end
    
    % Re-index D into a vector and sort the gammas into a matrix V such that
    % D is orthonormal to the columns of V when the network is bistochastic

    D_vec = [zeros(n*(n-1)/2,1);Const];
    ind=0;
    V = [zeros(n*(n-1)/2,2*n); -ones(1,2*n)]; % V dot D_vec will be zero for each column of V
    for i=1:n-1
        for j=1+i:n
            ind = ind + 1;
            D_vec(ind) = D_init(i,j);
            for k = 1:n
                if k==i
                    V(ind,k) = gamma(k,j);
                    V(ind,k+n) = gamma(j,k);
                elseif k==j
                    V(ind,k) = gamma(k,i);
                    V(ind,k+n) = gamma(i,k);
                end
            end
        end
    end
    
    % Step 3
    
    % Now find orthonormal basis for S^orthogonal
    Q = GramSchmidt(V);
    Q(:,rank(V):end) = [];
    % Projector matrix onto S^orthgonal
    P_orth = Q*Q';

    % Step 4
    % P = I - P^orth is the projector matrix onto S
    P = eye(n*(n-1)/2+1)-P_orth;

    % Find up to n*(n-1)/2+1-rank(v) columns of P with only positive components
    F = zeros(n*(n-1)/2+1,n*(n-1)/2+1-rank(V));
    count = 0;
    for i = 1:n*(n-1)/2+1
        if all(P(:,i)>0)
            count = count + 1;
            F(:,count) = P(:,i);
        end
        if count == n*(n-1)/2+1-rank(V)
            break
        end
    end
    
    % If there are no positive columns of P, skip this network
    if count == 0
        continue
    end

    % Step 5
    % Increase D's along those positive vectors in some way

    alpha = 1*rand(n*(n-1)/2+1-rank(V),1);
    D_vec = D_vec + F*diag(alpha);

    % Put increased D's back into a matrix
    D = zeros(n);
    ind=0;
    for i=1:n-1
        for j=1+i:n
            ind = ind + 1;
            D(i,j) = D_vec(ind);
        end
    end
    D = D+D';
    lambda = gamma.^2;
    C_init = (1+lambda)./(2*gamma).*D_init;
    C_final = (1+lambda)./(2*gamma).*D;
    
    R_eff(trial,1) = double(calc_R_eff(C_init,lambda));
    R_eff(trial,2) = double(calc_R_eff(C_final,lambda));
    
    trial = trial + 1;
    
end
R_eff
