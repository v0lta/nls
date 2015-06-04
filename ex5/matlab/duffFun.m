function [ dx ] = duffFun( t,x )
%A function that returns the change of variables for the duffing
%osscillator for use with ode45.

dx = zeros(2,1);

%chaos unique n=2 pair
%B = 7.5;
%k = 0.05;
%k = 0.2;
%k = 0.185;

%chaos unique
%B = 11.5;
%k = 0.2
%k = 0.185;

%no chaos;
k = 1;
B = 5;

dx(1) = x(2);
dx(2) = -k*x(2) - x(1)^3 + B * cos(t);


end

