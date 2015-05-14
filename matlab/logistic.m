
% dn = rN(1 - N/K)

% growth rate 
r = 1;

% carrying capacity.
K = 4;

%cell number:
N = 0:0.1:6;
dN = zeros(1,length(N));

for i = 1:length(N)
   dN(i) = r*N(i)*(1 - N(i)/K);
end

figure(1)
plot(N,dN)
xlabel('# cells (N)')
ylabel('population growth rate (dN)')
grid on

%solve with forward euler and plot different trajectories:
for N0 = 0:9;
    step = 0.01;
    tend = 10;

    N = ones(1,floor(tend/step));
    N(1) = N0;

    for i = 1:floor(tend/step)
       
        N(i+1) = N(i) + step*(r*N(i)*(1 - N(i)/K));
    
    end

    figure(2)
    plot(0:step:tend,N)
    xlabel('time (t)')
    ylabel('# cells (N)')
    hold on
end