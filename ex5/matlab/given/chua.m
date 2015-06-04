function [ dydx ] = chua( t,y )
%A function that returns the change of variables for chua's circuit for use
%with ode1.

n=3;
alpha = 9;
beta = 14.286;
a = -(1/7);
b = 2/7;

if (y(1,1) > 1) 
    h = b*y(1,1)  + a - b;
    dh = b;
elseif (y(1,1) < (-1))
    h = b*y(1,1) - a + b;
    dh = b;
else
    h = a*y(1,1);
    dh = a;
end

%   Integrate phase point (y(1,1)=x, y(2,1)=y, y(3,1)=z)
dydx(1,1)= alpha * ( y(2,1) - h);
dydx(2,1)= y(1,1) - y(2,1) + y(3,1);
dydx(3,1)= - beta * y(2,1);

%   Integrate tangent vectors with Jacobian (first component myin, second
%   component myin+1, third component myin+2, for myin =4,7,10);
%   for Lorenz the Jacobian is   |-dh*alpha  alpha   0|
%                                |   1       -1      1|
%                                |   0       -beta   0|

myin=n+1:n:n*(n+1);
dydx(myin,1)  =-dh*alpha*y(myin,1)+alpha*y(myin+1,1);
dydx(myin+1,1)=         +y(myin,1)      -y(myin+1,1)     +y(myin+2,1);
dydx(myin+2,1)=                    -beta*y(myin+1,1);
end

