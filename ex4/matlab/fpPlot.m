

x = 0:0.1:2;

a = 0.4;
b = 0.3;
c = 1.5;
d = 0;

fun1 = ((x-a).*(1-x))/b;
fun2 = d./(x-c);

plot(x,fun1)
hold on;
plot(x,fun2)