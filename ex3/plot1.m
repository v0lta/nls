
%%----------------------- find the zeros ----------------------------------

r = 0;
x = -5:0.1:5;
y =  0.5.*x.^3 - r.*x;

plot(x,y);
hold on;
h = -2;
plot(x,h);

coeff = [-0.5 0 r h];
rootVec = roots(coeff)





%------------------------plotting r with parameter h-----------------------


