function out = bazykin2
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
function dydt = fun_eval(t,kmrgd,a)
d=0;
e=0.01;
t=(kmrgd(1)*kmrgd(2))/(1+a*kmrgd(1));
dydt=[kmrgd(1)-t-d*kmrgd(1)^2;
-kmrgd(2)+t-e*kmrgd(2)^2;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(bazykin2);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,a)

% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,a)

% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,a)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,a)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,a)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,a)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,a)
