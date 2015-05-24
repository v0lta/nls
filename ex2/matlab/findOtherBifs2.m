
%% Entering parameters....
m = 1;
rho = 1;
r = 1;
k = 100;
a = 1;

A1 = 0.04695;
A3 = 8.932*10^(-4);
A5 = 1.015*10^(-5);
A7 = 2.955*10^(-8);

%% Plot:
V = 53.2481;

dy = 0:0.1:1000;

y = - dy/k + 0.5*V^2/k * (A1 .* dy/V - A3 .* (dy./V).^3 + A5.*(dy./V).^5 ...
        -A7.*(dy./V).^7);
    
plot(dy,y)
xlabel('dy')
ylabel('y')
