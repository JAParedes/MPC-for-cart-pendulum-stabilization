%Program implements Switching NMPC for extremely constrained
%wind-up maneuvers
%Scenario: Wind-up Maneuver with Small Position and Input Constraints
close all
clear all
clc

import casadi.*;

N = 30; %Prediction Horizon
nx = 4; %Number of states
nu = 1; %Number of inputs
%Parameters
ts = 0.1; %Sampling time
x0 = SX.sym('x0',nx);   %State Estimate
q = SX.sym('q',4,4);    %State weight matrix
r = SX.sym('r',1,1);    %Input weight matrix
qf = SX.sym('qf',4,4);  %Terminal state weight matrix
xTarget = SX.sym('xTarget',nx);   %Target Reference
xub = SX.sym('xub',nx,1);   %state upper bound
xlb = SX.sym('xlb',nx,1);   %state lower bound
uub = SX.sym('uub',nu,1);   %control upper bound
ulb = SX.sym('ulb',nu,1);   %control lower bound
gamma = SX.sym('gamma',1,1);   %penalty weight

%System parameters
m = 0.2;
M = 1;
L1 = 0.2;
J1 = m*L1*L1/3;
gm = 9.81;

%Creating stage functions
xk = SX.sym('xk',nx);   %State
uk = SX.sym('uk',nu);   %Input
sys =cartpendCas(0,xk,uk,m,M,L1,J1); %System dynamics equation
f_x = jacobian(sys,xk); %A matrix equation (linearized sys by x)
f_u = jacobian(sys,uk); %B matrix equation (linearized sys by u)
%*********************************************************
xkp1 = xk + ts*cartpendCas(0,xk,uk,m,M,L1,J1); %Discrete time dynamics obtained through
                                 %explicit Euler approximation
%*********************************************************

%Casadi function for discrete dynamics
fd = Function('f',{xk,uk},{xkp1},{'xk','uk'},'xkp1');
%Casadi function for obtaining linearized A matrix
f_x = Function('f_x',{xk,uk},{f_x},{'xbar','ubar'},{'f_x'});
%Casadi function for obtaining linearized B matrix
f_u = Function('f_u',{xk,uk},{f_u},{'xbar','ubar'},{'f_u'});
%Stage cost Casadi functions
ell = 0.5.*q(1,1).*(xk(1)-xTarget(1)).^2 + 0.5.*q(2,2).*(xk(2)-xTarget(2)).^2 + m.*gm.*q(3,3).*(cos(xTarget(3))-cos(xk(3))) + 0.5.*J1.*q(4,4).*(xk(4)-xTarget(4)).^2 + ...
    1/2*uk'*r*uk;
l = Function('l',{xk,uk,q,r,xTarget},{ell},{'xk','uk','q','r','xTarget'},{'l'});

%Creating primal optimization variables within finite horizon
x = SX.sym('x',nx,N);
u = SX.sym('u',nu,N);
s = SX.sym('s',N);
%Initializing J, h and g
h = [];
g = [];
J = SX.zeros(1);

%*********************************************************
%The following loop builds the J scalar and
%the h (inequalities) and g (equalities) vectors
%using z variables
for i = 1:N
   if i == 1
       xim1 = x0;
   else
       xim1 = x(:,i-1);
   end
   xi = x(:,i);
   uim1 = u(:,i);
   
   %cost
   J = J + l(xim1,uim1,q,r,xTarget) + s(i)*gamma;
   %inequality constraint
   h = [h;xi-s(i)-xub;-xi-s(i)+xlb;uim1-uub;-uim1+ulb;-s(i)];
   %equality constraint
   g = [g;xi-fd(xim1,uim1)];
   
end
%Terminal cost
J = J + 0.5.*qf(1,1).*(xi(1)-xTarget(1)).^2 + 0.5.*qf(2,2).*(xi(2)-xTarget(2)).^2 + m.*gm.*qf(3,3).*(cos(xTarget(3))-cos(xi(3))) + 0.5.*J1.*qf(4,4).*(xi(4)-xTarget(4)).^2;
%**********************************************************
%Reshaping every matrix into a column
x = x(:);
u = u(:);
z = [u;x;s];

nz = length(z);
nv = length(h);
nl = length(g);

%Defining dual variables and Lagrangian
clear l
l = SX.sym('l',nl,1);   %Equality duals
v = SX.sym('v',nv,1);   %Inequality duals
L = J + l'*g + v'*h;    %Lagrangian

%**************************************************************
%On the following lines, the gradients and Hessians for functions
% in question 1d) are obtained through Casadi
%Differentiating functions
Lz = gradient(L,z); %Lagrangian gradient (vector)
Jzz = hessian(J,z); %Cost function hessian (matrix)
Jz = gradient(J,z); %Cost function gradient (vector)
gz = jacobian(g,z); %Equality constraints jacobian (matrix)
hz = jacobian(h,z); %Inequality constraints jacobian (matrix)
%**************************************************************

%*************************************************************
%After that, all of the equations from 2d) will be stored in 
%struct nlp for later use
%Laplacian (L) gradient with respect to z
nlp.laplace_grad = ...
    Function('dL',{x0,x,u,s,l,v,xub,xlb,uub,ulb,...
    q,r,qf,xTarget,gamma},{Lz},...
    {'x0','x','u','s','l','v','xub','xlb','uub','ulb',...
    'q','r','qf','xTarget','gamma'},{'Lz'});
%Cost function (J) Hessian with respect to z
nlp.cost_hess = ...
    Function('HJ',{q,r,qf,x,xTarget},{Jzz},{'q','r','qf','x','xTarget'},{'Jzz'});
%Cost function (J) Gradient with respect to z
nlp.cost_grad = ...
    Function('dJ',{x,u,q,r,qf,xTarget,gamma},{Jz},...
    {'x','u','q','r','qf','xTarget','gamma'},{'Jz'});
%Inequality constraint function (h) gradient with respect
%to z (results in matrix)
nlp.in_grad = ...
    Function('HG',{},{hz},{},{'hz'});
%Inequality constraint function (h) with respect
%to z (results in vector)
nlp.in_reg = ...
    Function('dG',{x,u,s,xub,xlb,uub,ulb},{h},...
    {'x','u','s','xub','xlb','uub','ulb'},{'h'});
%Equality constraint function (g) gradient with respect
%to z (results in matrix)
nlp.eq_grad = ...
    Function('Hh',{x,u,x0},{gz},...
    {'x','u','x0'},{'gz'});
%Inequality constraint function (h) with respect
%to z (results in vector)
nlp.eq_reg = ...
    Function('dh',{x0,x,u},{g},...
    {'x0','x','u'},{'g'});
%*************************************************************
%Maximum number of iterations in OCP
args.max_sqp_iters = 20;
%Tolerance for ending OCP program
args.tol = 1e-4;

xtrg = zeros(4,1); %Wind up maneuver target conditions

% Initial Conditions
x0 = [0;0;pi;0]; %Wind up maneuver

% state constraints

args.xub = [1, 5, 20, 5]; %Wind up hard constraints
args.xlb = -args.xub;
% control constraints
umax = 10;
args.uub = umax*ones(1,1);
args.ulb = -args.uub;

args.N = N;         %Finite horizon
args.nlp = nlp;     %Struct with OCP required functions

args.Q = diag([30,5,60,-0.005]); %First state cost matrix
args.Q2 = diag([20,10,60,10]);  %Second state cost matrix

args.R = 10;        %New input cost matrix
%Obtaining Terminal state cost matrix by solving DARE using linearized
%nonlinear system
[~,Qf,~] = dlqr(full(f_x(0,0)),full(f_u(0,0)),args.Q2,args.R);
%args.Qf = diag([30,0,60,-0.005]); %Terminal state cost matrix
args.Qf = diag([40,0,60,-0.005]); %First Terminal state cost matrix
args.Qf2 = Qf; %Second Terminal state cost matrix
args.gamma = 10; %Weight for s

% initial optimizer guess for u, x, s (z), l, v (duals)
ukm1 = zeros(1*N,1);
xkm1 = zeros(4*N,1);
skm1 = zeros(N,1);
lkm1 = zeros(4*N,1);
vkm1 = zeros(11*N,1);
tfinal = 10; %Total simulation time
nsim = ceil(tfinal/ts) + 1; %Number of iterations for simulation given time step

%Matrices for storing time, state and input data from simulation in case 1
X = zeros(4,nsim);
t1 = zeros(nsim,1);
U = zeros(1,nsim);
X(:,1) = x0;
x = x0;
u_nmpc = 0;
t1(1) = 0;

for i = 1:nsim-1

    %Switching NMPC iterative solver
    [u_nmpc,ukm1,xkm1,skm1,lkm1,vkm1,res1(i)] = ...
    switching_nmpc(x,xtrg,ukm1,xkm1,skm1,lkm1,vkm1,args);
    
    %If input obtained is unfeasible, then reset duals and
    %optimization varuables and move cart towards center
    if u_nmpc > umax && x(1)<=0
        u_nmpc = umax;
        xkm1 = zeros(4*N,1);
        ukm1 = zeros(1*N,1);
        skm1 = zeros(N,1);
        lkm1 = zeros(4*N,1);
        vkm1 = zeros(11*N,1);
        ukm1(1) = umax;

    elseif u_nmpc > umax && x(1)>0
        u_nmpc = -umax;
        xkm1 = zeros(4*N,1);
        ukm1 = zeros(1*N,1);
        skm1 = zeros(N,1);
        lkm1 = zeros(4*N,1);
        vkm1 = zeros(11*N,1);
        ukm1(1) = -umax;
    
    elseif u_nmpc < -umax
        u_nmpc = -umax;
        xkm1 = zeros(4*N,1);
        ukm1 = zeros(1*N,1);
        skm1 = zeros(N,1);
        lkm1 = zeros(4*N,1);
        vkm1 = zeros(11*N,1);
        ukm1(1) = -umax;
    end
    
    [tx,xm] = ode45(@(t,x) cartpendCas(t,x,u_nmpc,m,M,L1,J1),0:ts/10:ts,x);
    
    x = xm(end,:)';

    X(:,i+1) = x;
    U(:,i) = u_nmpc;
    t1(i) = (i-1)*ts;

    percentage_done = i/nsim
end

%X(3,:) = wrapToPi(X(3,:));

U(:,nsim) = U(:,nsim-1);
res1(nsim) = 0;
t1(nsim) = (nsim-1)*ts;

subplot(3,1,1)
plot(t1,X(1,:), 'LineWidth',2)
title('Nonlinear Switching MPC')
xlabel('t (sec)')
ylabel('Cart Position (m)')
legend('x_1(t)')

subplot(3,1,2)

plot(t1,X(3,:), 'LineWidth',2)
xlabel('t (sec)')
ylabel('Pendulum Angle (rad)')
legend('x_3(t)')

subplot(3,1,3)
plot(t1,U,'LineWidth',2)
xlabel('t (sec)')
ylabel('Horizontal Force (N)')
legend('u(t)')

%% Saving for animation (uncomment)

% Ts = ts;
% pos = X(1,:);
% vel = X(2,:);
% ang = X(3,:);
% angvel = X(4,:);
% inp = U;
% tt = t1.';
% 
% save('Inverted_Pendulum_on_a_Cart_Swing_Up_Switching_MPC.mat','Ts','pos','vel','ang','angvel','inp','tt')

%End Program