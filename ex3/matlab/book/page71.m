
%% ----------------------- find the zeros ----------------------------------

maxR = 30;
h = 5;

figure(1)
x = -5:0.1:5;
for r = linspace(-maxR,maxR,10)
    plot(x,y);
    hold on;
    y = -x.^3 + r.*x;
end
plot(x,h,'r.');
hold off;
axis([min(x),max(x),-100,100])

coeff = [-1 0 r h];
rootVec = roots(coeff)

%% ------------------------plotting r with parameter h-----------------------
figure(2)

for r = -maxR:0.1:maxR
    coeff = [-1 0 r h];
    rootVec = roots(coeff);
    for i = 1:length(rootVec)
        if isreal(rootVec(i))
            plot(r,rootVec(i),'b.')
        end
    end
    hold on;
end
hold off;