function out = fitzhugh
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
function dydt = fun_eval(t,kmrgd,I)
dydt=[10*(kmrgd(1)-kmrgd(1)^3/3-kmrgd(2)+I);
0.8*(-kmrgd(2)+1.25*kmrgd(1)+1.5);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(fitzhugh);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,I)
jac=[[10-10*kmrgd(1)^2,-10];[1,-4/5]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,I)
jacp=[[10];[0]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,I)
hess1=[[-20*kmrgd(1),0];[0,0]];
hess2=[[0,0];[0,0]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,I)
hessp1=[[0,0];[0,0]];
hessp(:,:,1) =hessp1;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,I)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,I)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,I)
