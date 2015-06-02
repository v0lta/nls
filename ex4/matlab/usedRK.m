
figure(1)




[T,Y2] = ode45(@predPreySys,[0 400],[1 0.3]);
plot(Y2(:,1),Y2(:,2))
hold on;

[T,Y3] = ode45(@predPreySys,[0 400],[0.67 0.3]);
plot(Y3(:,1),Y3(:,2))
hold on;

[T,Y1] = ode45(@predPreySys,[0 400],[1 0.4]);
plot(Y1(:,1),Y1(:,2))
hold on;

plot(Y2(1,1),Y2(1,2),'*')
plot(Y3(1,1),Y3(1,2),'*')
plot(Y1(1,1),Y1(1,2),'*')

hold off;


figure(2)
[T2,Y2] = ode45(@predPreySys,[0 400],[1 0.3]);
plot(T2,Y2(:,1))
hold on;

[T3,Y3] = ode45(@predPreySys,[0 400],[0.67 0.3]);
plot(T3,Y3(:,1))

[T1,Y1] = ode45(@predPreySys,[0 400],[1 0.4]);
plot(T1,Y1(:,1))
hold on;


plot(T2,Y2(:,2))
plot(T3,Y3(:,2))
plot(T1,Y1(:,2))
hold off