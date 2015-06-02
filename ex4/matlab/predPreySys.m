function [ dx ] = predPreySys( t,x)
%ode45 compatible function for the fourth ex(1)ercise.

dx = zeros(2,1);

a = 0.4;
b = 0.3;
c = 0.75;
d = 0;

dx(1) = x(1)*(x(1) - a)*(1 - x(1)) - b*x(1)*x(2);
dx(2) = x(1)*x(2) - c*x(2) - d;

end

