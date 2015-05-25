function [ dy ] = bridge3( t,y)
%This function represents the system of first order equations for the 
%second nls exercise in input and output vectors contain the derivatives,
%from the zeroth' to the second derivative in increasing order. 

dy = zeros(3,1);

m = 1;
rho = 1;
r = 1;
k = 100;
a = 1;

A1 = 0.04695;
A3 = 8.932*10^(-4);
A5 = 1.015*10^(-5);
A7 = 2.955*10^(-8);

V = y(3);
%dot(z) = inVec(3); z = inVec(2); y = inVec(1);

dy(1) = y(2);
dy(2) = -k/m * y(1) - r/m * y(2) + (0.5/m)*(rho*V^2*a) ...
            * (A1 * (y(2)/V)^1 - A3 * (y(2)/V)^3 + A5 * (y(2)/V)^5 ...
            - A7 * (y(2)/V)^7);

end

