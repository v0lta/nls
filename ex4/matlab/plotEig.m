

a = 0.4;
b = 0.3;

c = 0:0.05:1.5;

figure(1)
for i = 1:length(c)
    J = [c(i)*(1 + a - 2*c(i))  -b*c(i); 
         (c(i)-a)*(1-c(i))/b    0]; 
    eigs = eig(J);
    plot(real(eigs),imag(eigs),'*');
    hold on;
end
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

figure(2)
for i = 1:length(c)
    J = [c(i)*(1 + a - 2*c(i))  -b*c(i); 
         (c(i)-a)*(1-c(i))/b    0]; 
    eigs = eig(J);
    if c(i) < 0.5
        plot(real(eigs),imag(eigs),'b*');
        hold on;
    end
    if ((c(i) > 0.5) && ( c(i) < 1))
        plot(real(eigs),imag(eigs),'r*');
        hold on;
    end
    if ( c(i) > 1)
        plot(real(eigs),imag(eigs),'g*');
        hold on;
    end
end
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')



figure(3)
eigs = zeros(2,length(c));
for i = 1:length(c)
    J = [c(i)*(1 + a - 2*c(i))  -b*c(i); 
         (c(i)-a)*(1-c(i))/b    0]; 
    eigs(:,i) = eig(J);
    hold on;
end
plot3(c,real(eigs),imag(eigs))
xlabel('c')
ylabel('Re(\lambda)')
zlabel('Im(\lambda)')


figure(4)
subplot(1,2,1)
plot(c,real(eigs))
xlabel('c')
ylabel('Re(\lambda)')
subplot(1,2,2)
plot(c,imag(eigs))
xlabel('c')
ylabel('Im(\lambda)')


