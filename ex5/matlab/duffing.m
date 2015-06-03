

[T,X] = ode45(@duffFun,[0 200],[0 0]);

plot(X(:,1),X(:,2))