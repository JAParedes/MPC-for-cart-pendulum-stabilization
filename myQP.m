function [U, lam] = myQP(H, q, A, b, lam0)
    %Simple dual projected optimization algorithm for solving Quadratic
    %Programming (QP) problems.
    %Max number of iterations
    maxIt = 100;

    %Getting matrices Hd and qd
    Hinv = inv(H);
    A_Hinv = A*Hinv;
    Hd = A_Hinv*A';
    qd = A_Hinv*q + b;

    %Obtaining the Lyapunov constant
    L = norm(Hd,inf);
    %Setting the initial lambda guess
    lam  = lam0;

    %Iterating gradient descent until convergence (maxIt)
    for i = 1:maxIt
        lam = max(0,lam - (1/L)*(Hd*lam + qd));
    end
    %Obtaining minimizing values
    U = -Hinv*(q+A'*lam);
end