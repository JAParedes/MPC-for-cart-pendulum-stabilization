function xdot = cartpend(t, x, u, m, M, L, J)

%Cart-Pendulum System Dynamcs
%Parameters
g = 9.80655;
a = (J + m*L^2);
b = (m*L)^2;

%Cart velocity
xdot(1) = x(2);
%Pendulum angular velocity
xdot(3) = x(4);
%Cart acceleration
xdot(2) = (a*m*L*sin(x(3))*(x(4))^2 - b*g*sin(x(3))*cos(x(3)) + a*u)/((m+M)*a - b*cos(x(3)));
%Pendulum angular acceleration
xdot(4) = ((m+M)*m*g*L*sin(x(3)) - b*sin(x(3))*cos(x(3))*(x(4))^2 - m*L*cos(x(3))*u)/((m+M)*a - b*cos(x(3)));

xdot = xdot';

end

