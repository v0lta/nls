function out = lorentz1
%
% Standard ode file of Lorentz
%
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];%@hessians;
out{6} = [];%@hessiansp;
out{7} = [];%@der3
out{8} = [];%@der4
out{9} = [];%@der5

% --------------------------------------------------------------------------
function dydt = fun_eval(t,x,sigma,r,b)
dydt = [sigma*(-x(1)+x(2)); r*x(1)-x(2)-x(1)*x(3); -b*x(3)+x(1)*x(2)];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
tspan = [0; 6.19216933131963970674];
y0 = [0;50;600];
handles = feval(@lorentz);
options = odeset('Jacobian',handles(3),'JacobianP',handles(4));
% --------------------------------------------------------------------------
function jac = jacobian(t,x,sigma,r,b)
jac = [[-sigma sigma 0];[r-x(3) -1 -x(1)];[x(2) x(1) -b]];

% --------------------------------------------------------------------------
function jacp = jacobianp(t,x,sigma,r,b)
jacp = [[-x(1)+x(2) 0 0];[0 x(1) 0];[0 0 -x(3)]];
