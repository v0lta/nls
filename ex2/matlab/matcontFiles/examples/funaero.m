function out = funaero
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,f)
A1 = 0.04695;
A3 = 8.932e-4;
A5 = 1.015e-5;
A7 = 2.955e-8;
Vc = 2/A1;
V  = Vc*f;

alfa = kmrgd(2)/f;
dydt=[kmrgd(2);
    -kmrgd(2) - 100*kmrgd(1) + 0.5*(V^2/Vc)*(A1*alfa - A3*alfa^3 + A5*alfa^5 - A7*alfa^7);];
% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(aeronlin);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,f)

% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,f)

% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,f)

% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,f)

%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,f)

%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,f)

%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,f)

