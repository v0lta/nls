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
function dydt = fun_eval(t,x,r,h)

dydt = -0.5*x^3 + r*x + h;

% -------------------------------------------------------------------------
function [tspan,y0,options] = init
y0=[0,0];

options = odeset;
handles = feval(myfuncont);
tspan = [0 10];

% -------------------------------------------------------------------------
function jac = jacobian(t,x,r,h)
