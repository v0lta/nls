
%%Entering parameters....
m = 1;
rho = 1;
r = 1;
k = 100;
a = 1;

A1 = 0.04695;
%A3 = 8.932*10^(-4);
%A5 = 1.015*10^(-5);
%A7 = 2.955*10^(-8);

%% Computing Vc:
Vc = (2*r)/(rho*a*A1)


V = -Vc:0.1:2*Vc;
for i= 1:length(V) 
    tau1(i) = 0.5*(rho*V(i)*a*A1)/(m) - r/m;
    delta1(i) = (k/m);
end 

plot(V,tau1)
hold on;
plot(Vc,0,'*')
grid on;
xlabel('speed V')
ylabel('trace')
%plot(V,delta1)    

