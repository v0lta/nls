function out = pitchfork
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,r)
dydt=[r*kmrgd(1)-kmrgd(1)^3;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(pitchfork);
y0=[0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,r)
jac=[ r - 3*kmrgd(1)^2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,r)
jacp=[ kmrgd(1) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,r)
hess1=[ -6*kmrgd(1) ];
hess(:,:,1) =hess1;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,r)
hessp1=[ 1 ];
hessp(:,:,1) =hessp1;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,r)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,r)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,r)

