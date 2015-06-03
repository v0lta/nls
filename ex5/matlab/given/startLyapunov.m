
%
%x = [0 1 2500]'
x = [0 1 20]';
%x = [0.5 0.5 0.5]';
%x = [0.1 0.1 0.1]';
%x = [1 1 1]'
st = 0.01;
kkmax = 1000;
lyap = lyapunov(@rhs_lorenz,st,kkmax,x)

[T,Y] = ode45(@rhs_lorenzClone,[0 st*kkmax],x);
figure(2)
plot3(Y(:,1),Y(:,2),Y(:,3))