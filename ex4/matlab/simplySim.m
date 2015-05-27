
x0 = 0:0.1:1.5;

figure(1)
for i = 1:length(x0)
    [T,Y] = ode45(@simpleSys,[0 10],x0(i));
    plot(T,Y)
    hold on;    
end
xlabel('time [s]')
ylabel('x')
hold off;