
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

t = [0 20];
%bc = [0.5 150 Vc*1.0];
%bc = [0.5 0 Vc*1.1];
bc = [0.5 0 Vc*0.9];


%% Non-Linear simulation 1. ------------------------------------------
t = [0 50];
f = 0.9

bc = [0.5 1500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

figure(1)
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')
hold off;

%% Non-Linear simulation 2. -------------------------------------------
f = 1

bc = [0 5 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

figure(2)
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')

bc = [0 10 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')

bc = [0 18 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')


%% Non-Linear simulation 3. -------------------------------------------
t = [0 75];

f = 1.241
bc = [0 50 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

figure(3)
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')

bc = [0.5 1500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')

hold off;


%% Non-Linear simulation 4. -------------------------------------------
t = [0 50];
f = 1.25

bc = [0 50 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

figure(4)
plot(Y(:,1),Y(:,2),':',bc(1),bc(2),'*')
hold on;

bc = [0.5 500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);
plot(Y(:,1),Y(:,2),':')
[T,Y] = ode45(@bridge3back,t,bc);
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')


bc = [0.5 1500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')
hold off;


%% Non-Linear simulation 5. -------------------------------------------
f = 1.5

bc = [0 50 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

figure(5)
plot(Y(:,1),Y(:,2),':',bc(1),bc(2),'*')
hold on;

bc = [0.5 500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);
plot(Y(:,1),Y(:,2),':')
[T,Y] = ode45(@bridge3back,t,bc);
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')


bc = [0.5 1500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')
hold off;

%% Non-Linear simulation 6. -------------------------------------------
f = 1.9

bc = [0 50 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

figure(6)
plot(Y(:,1),Y(:,2),':',bc(1),bc(2),'*')
hold on;

bc = [0.5 500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);
plot(Y(:,1),Y(:,2),':')
[T,Y] = ode45(@bridge3back,t,bc);
plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')


bc = [0.5 1500 Vc*f];
[T,Y] = ode45(@bridge3,t,bc);

plot(Y(:,1),Y(:,2),':')
hold on;
plot(bc(1),bc(2),'*')
hold off;


hold off;




