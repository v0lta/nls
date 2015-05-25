
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

t = [1 3000];
bc = [0.5 0];
%% Non-Linear simulation.
[T,Y] = ode45(@bridge,t,bc);

figure(1)
plot(1/10 * (1 + 0.1*T/10),Y(:,1))
xlabel('V/Vc')
ylabel('y')


%t = [1 100];
%bc = [0.5 0];
%% Non-Linear simulation.
%[T,Y] = ode45(@bridge,t,bc);

%figure(1)
%plot(T,Y(:,1))
%ylabel('y')

%figure(2)
%plot(Y(:,1),Y(:,2))


%a limit cycle forms = hopf bifurcation.