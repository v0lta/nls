function [ dx ] = Lorenz( t,x )
%A lorenz-strange-attractor function for simulation with ode45.

dx = zeros(3,1);

a = 10;
b = 28;
c = 8/3;

dx(1) = a * (x(2) - x(1));
dx(2) = x(1) * (b - x(3)) - x(2);
dx(3) = x(1)*x(2) - c * x(3);


end

