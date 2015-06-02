

x = 0:0.001:2;

a = 0.4;
b = 0.3;
c = 0.2;
d = -0.01;

fun1 = ((x-a).*(1-x))/b;
fun2 = d./(x-c);

plot(x,fun1)
hold on;
plot(x,fun2)
hold off;