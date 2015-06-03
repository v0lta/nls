function dydx=rhs_lorenzClone(t,y)
b = 8/3; sigma = 10; r = 28;
n=3;
%   Integrate phase point (y(1)=x, y(2)=y, y(3)=z)
dydx(1,1)=sigma*(y(2,1)-y(1,1));
dydx(2,1)=r*y(1,1)-y(1,1)*y(3,1)-y(2,1);
dydx(3,1)=y(1,1)*y(2,1)-b*y(3,1);
