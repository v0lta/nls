function out = myfuncont
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% -------------------------------------------------------------------------
function dxdt = fun_eval(t,x,c,d)
a = 0.4;
b = 0.3;
%c = 0.75;
%d = 0;

dxdt = [  x(1)*(x(1) - a)*(1 - x(1)) - b*x(1)*x(2);
          x(1)*x(2) - c*x(2) - d;];
% -------------------------------------------------------------------------
function [tspan,y0,options] = init
x0=[0,0];

options = odeset;
handles = feval(myfuncont);
tspan = [0 50];

% -------------------------------------------------------------------------
function jac = jacobian(t,x,r,h)
