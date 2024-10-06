function xdot = cartpendCas(t, x, u, m, M, L, J)

%Inverted Pendulum on a Cart System Dynamcs compatible with Casadi programming scheme
%Parameters
g = 9.80655;
a = (J + m*L^2);
b = (m*L)^2;

xdot = [x(2);...
    (a*m*L*sin(x(3))*(x(4))^2 - b*g*sin(x(3))*cos(x(3)) + a*u)/((m+M)*a - b*cos(x(3)));...
    x(4);...
    ((m+M)*m*g*L*sin(x(3)) - b*sin(x(3))*cos(x(3))*(x(4))^2 - m*L*cos(x(3))*u)/((m+M)*a - b*cos(x(3)))];

end