clear
%% globals
%small values
%setArray = 0.1; %<--has to be pos!
%setParam = 0.001;
%step = 0.0001

%big values
setArray = 30; %<--has to be pos!
setParam = 0;
step = 0.5;


%% ----------------------- find the zeros ---------------------------------

maxR = setArray;
h = setParam;

figure(1)
x = linspace(-setParam,setParam,100);
for r = linspace(-maxR,maxR,10)
    y = -x.^3 + r.*x;
    plot(x,y);
    hold on;
end
plot(x,h,'r.');
hold off;
axis([min(x),max(x),-100,100])

coeff = [-1 0 r h];
rootVec = roots(coeff)

%% ------------------------plotting r with const parameter h---------------
figure(2)
h = setParam;
for r = -maxR:step:maxR
    coeff = [-1 0 r h];
    rootVec = roots(coeff);
    for i = 1:length(rootVec)
        if isreal(rootVec(i))
            plot(r,rootVec(i),'b.')
        end
    end
    hold on;
end
xlabel('r')
ylabel('roots(u)')
hold off;


%% ------------------------plotting h with const parameter r --------------- 
maxh = setArray;
r = setParam;

figure(3)
for h = -maxh:step:maxh
    coeff = [-1 0 r h];
    rootVec = roots(coeff);
    for i = 1:length(rootVec)
        if isreal(rootVec(i))
            plot(h,rootVec(i),'b.')
        end
    end
    hold on;
end
xlabel('h')
ylabel('roots(u)')
hold off;


