function out = MorrisLecar(t,coordinates,flag,epsilon,k)
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
function dydt = fun_eval(t,kmrgd,epsilon,k)
minf=(1+tanh((kmrgd(1)+0.01)/0.15))/2;
winf=(1+tanh((kmrgd(1)-0.1)/0.145))/2;
tau=cosh((kmrgd(1)-0.1)/0.29);
dydt=[kmrgd(3)-0.5*(kmrgd(1)+0.5)-2*kmrgd(2)*(kmrgd(1)+0.7)-minf*(kmrgd(1)-1);
1.15*(winf-kmrgd(2))*tau;
epsilon*(k-kmrgd(1));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(MorrisLecar);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,epsilon,k)
jac=[[-1-2*kmrgd(2)-(10/3-10/3*tanh(20/3*kmrgd(1)+1/15)^2)*(kmrgd(1)-1)-1/2*tanh(20/3*kmrgd(1)+1/15),-2*kmrgd(1)-7/5,1];[(115/29-115/29*tanh(200/29*kmrgd(1)-20/29)^2)*cosh(100/29*kmrgd(1)-10/29)+100/29*(23/40+23/40*tanh(200/29*kmrgd(1)-20/29)-23/20*kmrgd(2))*sinh(100/29*kmrgd(1)-10/29),-23/20*cosh(100/29*kmrgd(1)-10/29),0];[-epsilon,0,0]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,epsilon,k)
jacp=[[0,0];[0,0];[k-kmrgd(1),epsilon]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,epsilon,k)
hess1=[[20/3*tanh(20/3*kmrgd(1)+1/15)*(20/3-20/3*tanh(20/3*kmrgd(1)+1/15)^2)*(kmrgd(1)-1)-20/3+20/3*tanh(20/3*kmrgd(1)+1/15)^2,-2,0];[-230/29*tanh(200/29*kmrgd(1)-20/29)*(200/29-200/29*tanh(200/29*kmrgd(1)-20/29)^2)*cosh(100/29*kmrgd(1)-10/29)+200/29*(115/29-115/29*tanh(200/29*kmrgd(1)-20/29)^2)*sinh(100/29*kmrgd(1)-10/29)+10000/841*(23/40+23/40*tanh(200/29*kmrgd(1)-20/29)-23/20*kmrgd(2))*cosh(100/29*kmrgd(1)-10/29),-115/29*sinh(100/29*kmrgd(1)-10/29),0];[0,0,0]];
hess2=[[-2,0,0];[-115/29*sinh(100/29*kmrgd(1)-10/29),0,0];[0,0,0]];
hess3=[[0,0,0];[0,0,0];[0,0,0]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,epsilon,k)
hessp1=[[0,0,0];[0,0,0];[-1,0,0]];
hessp2=[[0,0,0];[0,0,0];[0,0,0]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,epsilon,k)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,epsilon,k)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,epsilon,k)
