%Program compares regular LQR against LQR with the reference signal given by a reference governor.
%Scenario: Cart movement while keeping pendulum upwards

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
%Obtaining dlqr gain
[K,~,~] = dlqr(Ad,Bd,Q,R);

%Setting system constraints
xlim.max = [22 100 100 100];
xlim.min = -xlim.max;
%umax = 100; %Unconstrained
umax = 20;
ulim.max = umax;
ulim.min = -ulim.max;

%Setting prediction horizon
N = 60;

%Setting matrices for Reference Governor implementation
Anew = Ad-Bd*K;
Bnew = Bd*K;
Cnew = [1 0 0 0;...
    0 0 1 0];
Dnew = zeros(2,4);
H1 = eye(2);
ub = [102;0.4];
lb = -ub;
eps = 0.01;

%[Qmat,hvec] = getMats(Anew, Bnew, Cnew, Dnew, H1, lb, ub, N, eps);
%Get matrices for evaluating whether future states will be within
%approximated safe zone
[Qmat,hvec] = getMats(Anew, Bnew, Cnew, Dnew, H1, K, lb, ub, -umax, umax, N, eps);

t = 0:Ts:15;
h = Ts/10;

%For storing results from LQR w/ reference governor implementation
X = zeros(4,length(t));
U = zeros(1,length(t));

%For storing results from LQR implementation
Xk = zeros(4,length(t));
Uk = zeros(1,length(t));

%Reference vector
r = [100;0;0;0];
%Initial conditions
x = [0;0;0;0];
tm = 0;

X(:,1) = x;
Xk(:,1) = x;
u = 0;%Input
v = 0; %Governed reference
Kvar = 1; %Used to determine how similar v is to r

%Reference governor implementation
for i = 1:length(t)-1
   %Searching for Kvar value that prevents constraint violation
   Kvar = 1;
   while 1
      vn = v + Kvar*(r-v);
      if prod(Qmat*[vn;x]-hvec<=0)==1
         v = vn;
         break;
      else
         Kvar = Kvar - 0.05; 
      end
       
      if Kvar <=0
          break;
      end
      
   end
   
   %Implementing new reference
   u = K*(v-x);
   
   [tx,xm] = ode45(@(t,x) cartpend(t,x,u,m,M,L1,J), tm+h:h:tm+Ts,x);
   x = xm(end,:)';
   %Keeps angle bounded for the LQR
   if x(3)>pi
        x(3)= x(3) - 2*pi;
    elseif x(3)<=-pi
        x(3)= x(3) + 2*pi;
    end
   U(i) = u;
   X(:,i+1) = x;
end

x = [0;0;0;0];
tm = 0;

%LQR implementation
for i = 1:length(t)-1
   u = K*(r-x);
   
   [tx,xm] = ode45(@(t,x) cartpend(t,x,u,m,M,L1,J), tm+h:h:tm+Ts,x);
   x = xm(end,:)';
   %Keeps angle bounded for the LQR
   if x(3)>pi
        x(3)= x(3) - 2*pi;
    elseif x(3)<=-pi
        x(3)= x(3) + 2*pi;
    end
   Uk(i) = u;
   Xk(:,i+1) = x;
end

subplot(3,2,1)
plot(t,X(1,:), 'LineWidth',2)
title('LQR w/ Reference Governor')
xlabel('t (sec)')
ylabel('Cart position (m)')
legend('x_1(t)')

subplot(3,2,3)
plot(t,X(3,:), 'LineWidth',2)
xlabel('t (sec)')
ylabel('Pendulum Angle (rad)')
legend('x_3(t)')

subplot(3,2,5)
plot(t,U,'LineWidth',2)
xlabel('t (sec)')
ylabel('Horizontal Force (N)')
legend('u(t)')

subplot(3,2,2)
plot(t,Xk(1,:), 'LineWidth',2)
title('LQR')
xlabel('t (sec)')
legend('x_1(t)')

subplot(3,2,4)
plot(t,Xk(3,:), 'LineWidth',2)
xlabel('t (sec)')
legend('x_3(t)')

subplot(3,2,6)
plot(t,Uk,'LineWidth',2)
xlabel('t (sec)')
legend('u(t)')

%Program end

% function [Q,h] = getMats(A,B,C,D,H,lb,ub,N,eps)
% 
%     Q = [H*D H*C; -H*D -H*C];
%     h = [ub;-lb];
%     invA = inv(eye(size(A))-A);
%     for i = 1:N
%         Q = [Q; H*C*(eye(size(A))-A^i)*invA*B+H*D H*C*A^i;...
%             -H*C*(eye(size(A))-A^i)*invA*B+H*D -H*C*A^i];
%         h = [h;ub;-lb];
%     end
%     Q = [Q;H*C*invA*B+H*D zeros(size(H*C));...
%         -H*C*invA*B+H*D -zeros(size(H*C))];
%      h = [h;(1-eps)*ub;-(1-eps)*lb];
% end

function [Q,h] = getMats(A,B,C,D,H,K,lb,ub,ulb,uub,N,eps)
    %Obtains matrices to test whether the future states considered within
    %the prediction horizon will remain inside the specified safe zone
    %under the given modified reference v
    Q = [H*D H*C; -H*D -H*C];
    h = [ub;-lb];
    invA = inv(eye(size(A))-A);
    for i = 1:N
        Q = [Q; H*C*(eye(size(A))-A^i)*invA*B+H*D H*C*A^i;...
            -H*C*(eye(size(A))-A^i)*invA*B+H*D -H*C*A^i];
        h = [h;ub;-lb];
    end
    Q = [Q;H*C*invA*B+H*D zeros(size(H*C));...
        -H*C*invA*B+H*D -zeros(size(H*C));...
        K -K; -K K];
    h = [h;(1-eps)*ub;-(1-eps)*lb; uub; -ulb];
end