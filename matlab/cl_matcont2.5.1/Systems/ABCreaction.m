function out = ABCreaction(t,coordinates,flag,p1,p2,p3,p4,p5)
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
function dydt = fun_eval(t,kmrgd,p1,p2,p3,p4,p5)
dydt=[-kmrgd(1)+p1*(1-kmrgd(1))*exp(kmrgd(3));
-kmrgd(2)+p1*(1-kmrgd(1)-p5*kmrgd(2))*exp(kmrgd(3));
-kmrgd(3)-p3*kmrgd(3)+p1*p4*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(ABCreaction);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,p1,p2,p3,p4,p5)
jac=[[-1-p1*exp(kmrgd(3)),0,p1*(1-kmrgd(1))*exp(kmrgd(3))];[-p1*exp(kmrgd(3)),-1-p1*p5*exp(kmrgd(3)),p1*(1-kmrgd(1)-p5*kmrgd(2))*exp(kmrgd(3))];[-p1*p4*exp(kmrgd(3)),p1*p4*p2*p5*exp(kmrgd(3)),-1-p3+p1*p4*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3))]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,p1,p2,p3,p4,p5)
jacp=[[(1-kmrgd(1))*exp(kmrgd(3)),0,0,0,0];[(1-kmrgd(1)-p5*kmrgd(2))*exp(kmrgd(3)),0,0,0,-p1*kmrgd(2)*exp(kmrgd(3))];[p4*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3)),p1*p4*p5*kmrgd(2)*exp(kmrgd(3)),-kmrgd(3),p1*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3)),p1*p4*p2*kmrgd(2)*exp(kmrgd(3))]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,p1,p2,p3,p4,p5)
hess1=[[0,0,-p1*exp(kmrgd(3))];[0,0,-p1*exp(kmrgd(3))];[0,0,-p1*p4*exp(kmrgd(3))]];
hess2=[[0,0,0];[0,0,-p1*p5*exp(kmrgd(3))];[0,0,p1*p4*p2*p5*exp(kmrgd(3))]];
hess3=[[-p1*exp(kmrgd(3)),0,p1*(1-kmrgd(1))*exp(kmrgd(3))];[-p1*exp(kmrgd(3)),-p1*p5*exp(kmrgd(3)),p1*(1-kmrgd(1)-p5*kmrgd(2))*exp(kmrgd(3))];[-p1*p4*exp(kmrgd(3)),p1*p4*p2*p5*exp(kmrgd(3)),p1*p4*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3))]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,p1,p2,p3,p4,p5)
hessp1=[[-exp(kmrgd(3)),0,(1-kmrgd(1))*exp(kmrgd(3))];[-exp(kmrgd(3)),-p5*exp(kmrgd(3)),(1-kmrgd(1)-p5*kmrgd(2))*exp(kmrgd(3))];[-p4*exp(kmrgd(3)),p4*p2*p5*exp(kmrgd(3)),p4*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3))]];
hessp2=[[0,0,0];[0,0,0];[0,p1*p4*p5*exp(kmrgd(3)),p1*p4*p5*kmrgd(2)*exp(kmrgd(3))]];
hessp3=[[0,0,0];[0,0,0];[0,0,-1]];
hessp4=[[0,0,0];[0,0,0];[-p1*exp(kmrgd(3)),p1*p2*p5*exp(kmrgd(3)),p1*(1-kmrgd(1)+p2*p5*kmrgd(2))*exp(kmrgd(3))]];
hessp5=[[0,0,0];[0,-p1*exp(kmrgd(3)),-p1*kmrgd(2)*exp(kmrgd(3))];[0,p1*p4*p2*exp(kmrgd(3)),p1*p4*p2*kmrgd(2)*exp(kmrgd(3))]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,p1,p2,p3,p4,p5)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,p1,p2,p3,p4,p5)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,p1,p2,p3,p4,p5)
