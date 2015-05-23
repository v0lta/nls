
%%Entering parameters....
m = 1;
rho = 1;
r = 1;
k = 100;
a = 1;

A1 = 0.04695;
A3 = 8.932*10^(-4);
A5 = 1.015*10^(-5);
A7 = 2.955*10^(-8);

%% Computing Vc:
Vc = (2*r)/(rho*a*A1);

[T,Y] = ode45(@bridge,[0 10],[0.5 0]);

plot(T,Y(:,1))
hold on;
plot(T,Y(:,2))
