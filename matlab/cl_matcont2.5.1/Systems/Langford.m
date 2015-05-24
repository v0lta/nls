function out = Langford
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = [];
out{8} = [];
out{9} =[];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,a,b,d,lambda)
dydt=[(lambda - b)*kmrgd(1) - a*kmrgd(2) + kmrgd(1)*kmrgd(3) + d*kmrgd(1)*(1-kmrgd(3)^2);
a*kmrgd(1) + (lambda - b)*kmrgd(2) + kmrgd(2)*kmrgd(3) + d*kmrgd(2)*(1-kmrgd(3)^2);
lambda*kmrgd(3) - (kmrgd(1)^2 + kmrgd(2)^2 + kmrgd(3)^2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Langford);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,a,b,d,lambda)
jac=[[lambda-b+kmrgd(3)+d*(1-kmrgd(3)^2),-a,kmrgd(1)-2*d*kmrgd(1)*kmrgd(3)];[a,lambda-b+kmrgd(3)+d*(1-kmrgd(3)^2),kmrgd(2)-2*d*kmrgd(2)*kmrgd(3)];[-2*kmrgd(1),-2*kmrgd(2),lambda-2*kmrgd(3)]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,a,b,d,lambda)
jacp=[[-kmrgd(2),-kmrgd(1),kmrgd(1)*(1-kmrgd(3)^2),kmrgd(1)];[kmrgd(1),-kmrgd(2),kmrgd(2)*(1-kmrgd(3)^2),kmrgd(2)];[0,0,0,kmrgd(3)]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,a,b,d,lambda)
hess1=[[0,0,1-2*d*kmrgd(3)];[0,0,0];[-2,0,0]];
hess2=[[0,0,0];[0,0,1-2*d*kmrgd(3)];[0,-2,0]];
hess3=[[1-2*d*kmrgd(3),0,-2*d*kmrgd(1)];[0,1-2*d*kmrgd(3),-2*d*kmrgd(2)];[0,0,-2]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,a,b,d,lambda)
hessp1=[[0,-1,0];[1,0,0];[0,0,0]];
hessp2=[[-1,0,0];[0,-1,0];[0,0,0]];
hessp3=[[1-kmrgd(3)^2,0,-2*kmrgd(1)*kmrgd(3)];[0,1-kmrgd(3)^2,-2*kmrgd(2)*kmrgd(3)];[0,0,0]];
hessp4=[[1,0,0];[0,1,0];[0,0,1]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,a,b,d,lambda)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,a,b,d,lambda)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,a,b,d,lambda)
