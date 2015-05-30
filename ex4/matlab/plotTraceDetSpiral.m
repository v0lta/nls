
%given conditions.
a = 0.4;
b = 0.3;

%Learn more about the stability of the Predetor-Prey Fixed points by
%plotting trace determinant and spiral condition.

figure(1)
x1 = [0 0];
x3 = [1 0];
x4 = [a 0];

plot(x1(1),x1(2),'x')
hold on;
plot(x3(1),x3(2),'s')
plot(x4(1),x4(2),'d')

for c = 0:0.1:1.5
    x2 = [c ((c-1)*(1-c))/b];
    plot(x2(1),x2(2),'p')
end
title('fixed point posistion.')
xlabel('x')
ylabel('y')
hold off;

figure(2)
title('trace and determinant')
hold on;
grid on;
det = 0:0.01:1;
plot(det,sqrt(4*det),'k')
plot(det,-sqrt(4*det),'k')
plot(det,zeros(length(det),1),'k')
c = 0:0.05:1.5;
for i = 1:length(c)
    % x1:
    tau1(i) = -a - c(i);
    det1(i) = a*c(i);
   
    %x2:
    tau2(i) = c(i)*(1 + a - 2*c(i));
    det2(i) = c(i)*(c(i) - a)*(1 - c(i));
    
    %x3 :
    tau3(i) = -c(i) + a;
    det3(i) = (a-1)*(1-c(i));
    
    %x4
    tau4(i) = a^2 + 2*a -c(i);
    det4(i) = (-a^2 + a)*(a-c(i));
end
plot(det1,tau1,'-x')
plot(det2,tau2,'-p')
plot(det3,tau3,'-s')
plot(det4,tau4,'-d')    
xlabel('det')
ylabel('tau')
hold off;

figure(3)
title('trace')
hold on;
plot(c,tau1,'-x')
plot(c,tau2,'-p')
plot(c,tau3,'-s')
plot(c,tau4,'-d')
xlabel('c')
ylabel('tau')
grid on; hold off;

figure(4)
title('det')
hold on;
grid on;
plot(c,det1,'-x')
plot(c,det2,'-p')
plot(c,det3,'-s') %square
plot(c,det4,'-d') %diamond
xlabel('c')
ylabel('det')
grid on; hold off;

figure(5)
title('3d tau,det')
hold on;
grid on;
plot3(tau1,det1,c,'-x')
plot3(tau2,det2,c,'-p')
plot3(tau3,det3,c,'-s')
plot3(tau4,det4,c,'-d')    
xlabel('tau')
ylabel('det')
zlabel('c')
grid on; hold off;





