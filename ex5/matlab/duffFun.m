function [ dx ] = duffFun( t,x )
%A function that returns the change of variables for the duffing
%osscillator for use with ode45.

dx = zeros(2,1);

k = 0.2;
B = 12;

dx(1) = x(2);
dx(2) = -k*x(2) - x(1)^3 + B * cos(t);


end

