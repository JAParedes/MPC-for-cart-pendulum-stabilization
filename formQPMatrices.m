function [H, L, G, W, T, IMPC] = formQPMatrices(A, B, Q, R, P, xlim, ulim, N)
    %This function creates matrices to set MPC as a quadratic programming (QP) problem
    %Obtaining sizes of input and state
    nu = size(B,2);
    nx = size(A,2);
    %initial S matrix, diagonal full of B matrices
    S = kron(eye(N,N),B);
    %Initial M matrix, to be filled with factors of A
    M = [];
    for i = 1:N-1
        %This process shifts the diagonal terms
        %of the identity matrix by i
        S1 = zeros(N,N);
        S2 = eye(N-i,N-i);
        S1(i+1:end,1:N-i)=S2;
        %Used to shift the (A^i)*B diagonals
        S = S + kron(S1,(A^(i))*B);
        %A factors are continuously added
        M = [M;A^i];
    end
    M = [M;A^N]; %Last A factor is added
    %Building Qbar with Qs on its diagonal
    Qbar = kron(eye(N,N),Q);
    %Last value of Qbar is P
    Qbar((end+1-size(P,1)):end,(end+1-size(P,2)):end) = P;
    %Rbar with R matrices on its diagonals
    Rbar = kron(eye(N,N),R);
    
    %Obtaining equivalent H for QP problem
    H = 2*(S'*Qbar*S + Rbar);
    %Obtaining equivalent L for QP problem such that
    %q = L*x0
    L = 2*S'*Qbar'*M;
    
    %Obtaining equivalent G, W, T for inequality constraints
    G = [S;-S;eye(N*nu,N*nu);-eye(N*nu,N*nu)];
    W = [kron(ones(N,1),xlim.max');-kron(ones(N,1),xlim.min');kron(ones(N,1),ulim.max');-kron(ones(N,1),ulim.min')];
    T = [-M;M;zeros(2*N*nu,nx)];
    %u0 = IMPC*U
    IMPC = [eye(nu),zeros(nu,nu*(N-1))];

end