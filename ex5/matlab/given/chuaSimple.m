function [ dydx ] = chuaSimple( t,y )
%A function that returns the change of variables for chua's circuit for use
%with ode1.

%-----------interface-----------------------------------
% fromChua = chua(t,y);
% dydx(1,1)= fromChua(1,1);
% dydx(2,1)= fromChua(2,1);
% dydx(3,1)= fromChua(3,1);

n=3;
alpha = 9;
beta = 14.286;
a = -(1/7);
b = 2/7;

if (y(1,1) > 1) 
    h = b*y(1,1)  + a - b;
elseif (y(1,1) < (-1))
    h = b*y(1,1) - a + b;
else
    h = a*y(1,1);
end

%   Integrate phase point (y(1,1)=x, y(2,1)=y, y(3,1)=z)
dydx(1,1)= alpha * ( y(2,1) - h);
dydx(2,1)= y(1,1) - y(2,1) + y(3,1);
dydx(3,1)= - beta * y(2,1);


end

