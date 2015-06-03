

%compute values
tend = 20;
[T1, X1] = ode45(@Lorenz,[1 tend],[.2 0 0]);
[T2, X2] = ode45(@Lorenz,[1 tend],[.2 0 0.1]);
[T3, X3] = ode45(@Lorenz,[1 tend],[.2 0 0.2]);
[T4, X4] = ode45(@Lorenz,[1 tend],[.2 0 0.3]);

%polt them
%arrange some nice colors.
colorMat = lines;
color1 = colorMat(1,:);
color2 = colorMat(2,:);
color3 = colorMat(3,:);
color4 = colorMat(4,:);


plot3(X1(:,1),X1(:,2),X1(:,3),'color',color1)
hold on
plot3(X2(:,1),X2(:,2),X2(:,3),'color',color2)
plot3(X3(:,1),X3(:,2),X3(:,3),'color',color3)
plot3(X4(:,1),X4(:,2),X4(:,3),'color',color4)


plot3(X1(end,1),X1(end,2),X1(end,3),'o','color',color1)
plot3(X2(end,1),X2(end,2),X2(end,3),'o','color',color2)
plot3(X3(end,1),X3(end,2),X3(end,3),'o','color',color3)
plot3(X4(end,1),X4(end,2),X4(end,3),'o','color',color4)




%[T, X] = ode45(@Lorenz,[1 200],[-10 -10 25]);

%plot3(X(:,1),X(:,2),X(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
grid on;
hold off;

