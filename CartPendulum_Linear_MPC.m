%Program compares Linear MPC against LQR 
%Scenario: Stabilizing pendulum at upwards position from initial angle
%close to upwards position

close all
clear all
clc

%System sampling time
Ts = 0.01;

x = sym('x',[4 1],'real');
u = sym('u',[1 1],'real');

%System parameters
m = 0.2;
M = 1;
L1 = 0.2;
J = m*L1*L1/3;

%Obtaining linearized system
fun = cartpend(0,x,u,m,M,L1,J);

Ap = jacobian(fun,x);
Bp = jacobian(fun,u);

Ac = double(subs(Ap, [x;u], [zeros(size(x));zeros(size(u))]));
Bc = double(subs(Bp, [x;u], [zeros(size(x));zeros(size(u))]));

clear x u

Cc = eye(size(Ac));
Dc = zeros(size(Bc));

sysd = c2d(ss(Ac,Bc,Cc,Dc),Ts,'zoh');

%Discretized linearized state space matrices
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

%Establishing input and state cost matrices
Q = diag([1 1 5 5]);
R = 0.5;
%Obtaining dlqr gain and DARE solution for terminal state cost matrix
[K,P,~] = dlqr(Ad,Bd,Q,R);

%Setting system constraints
xlim.max = [20 300 300 300];
xlim.min = -xlim.max;
%umax = 100; %Unconstrained
umax = 30;
umin = -umax;
ulim.max = umax;
ulim.min = -ulim.max;

%Setting prediction horizon
N = 20;

%Obtaining quadratic programming matrices
[H, L, G, W, T, IMPC] = formQPMatrices(Ad, Bd, Q, R, P, xlim, ulim, N);
%Setting initial guess of lambda
lam0 = ones(size(G,1),1);
lam = lam0;

t = 0:Ts:10;
h = Ts/10;

%Storing results from Linear MPC
X = zeros(4,length(t));
U = zeros(1,length(t));

%Storing results for LQR
Xk = zeros(4,length(t));
Uk = zeros(1,length(t));

%Initial position
x0 = pi/2;

x = [0;0;x0;0];
tm = 0;

X(:,1) = x;
Xk(:,1) = x;

for i = 1:length(t)-1
   
   [du,lam] = myQP(H, L*x, G, W + T*x, lam);
    %Determine whether or not output matrix is empty
    if ~isempty(du)
        u = IMPC*du;
    else
        u = 0;
    end
    %Implementing input constraints
    if u > umax
        u = umax;
    elseif u < umin
       u = umin; 
    end
   %Simulating at each sampling time 
   [tx,xm] = ode45(@(t,x) cartpend(t,x,u,m,M,L1,J), tm+h:h:tm+Ts,x);
    x = xm(end,:)';
   %Keeps angle bounded for the linear MPC
   if x(3)>pi
        x(3)= x(3) - 2*pi;
    elseif x(3)<=-pi
        x(3)= x(3) + 2*pi;
    end
   U(i) = u;
   X(:,i+1) = x;
end

x = [0;0;x0;0];
tm = 0;

for i = 1:length(t)-1
   u = -K*x; %LQR linear feedback
   
   %Implementing input constraints
    if u > umax
        u = umax;
    elseif u < umin
       u = umin; 
    end
   %Simulating at each sampling time 
   [tx,xm] = ode45(@(t,x) cartpend(t,x,u,m,M,L1,J), tm+h:h:tm+Ts,x);
   x = xm(end,:)';
    %Keeps angle bounded for the linear feedback
   if x(3)>pi
        x(3)= x(3) - 2*pi;
    elseif x(3)<=-pi
        x(3)= x(3) + 2*pi;
    end
   Uk(i) = u;
   Xk(:,i+1) = x;
end

subplot(2,2,1)
plot(t,X, 'LineWidth',2)
title('Linear MPC')
xlabel('t (sec)')
legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'x_4(t)')

subplot(2,2,3)
plot(t,U,'LineWidth',2)
xlabel('t (sec)')
legend('u(t)')

subplot(2,2,2)
plot(t,Xk, 'LineWidth',2)
title('LQR')
xlabel('t (sec)')
legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'x_4(t)')

subplot(2,2,4)
plot(t,Uk,'LineWidth',2)
xlabel('t (sec)')
legend('u(t)')

%Program end