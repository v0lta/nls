function out = pdp2
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,alpha,beta)
dydt=[kmrgd(1)*(alpha-kmrgd(1)-6*kmrgd(2)-4*kmrgd(3));
kmrgd(2)*(beta-kmrgd(1)-kmrgd(2)-10*kmrgd(3));
-kmrgd(3)*(1-0.25*kmrgd(1)-4*kmrgd(2)+kmrgd(3));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(pdp2);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,alpha,beta)

% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,alpha,beta)

% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,alpha,beta)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,alpha,beta)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,alpha,beta)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,alpha,beta)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,alpha,beta)
