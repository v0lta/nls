function out = SteinmetzLarter
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= @test1;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
dydt=[-k1*kmrgd(1)*kmrgd(2)*kmrgd(3)-k3*kmrgd(1)*kmrgd(2)*kmrgd(4)+k7-km7*kmrgd(1);
-k1*kmrgd(1)*kmrgd(2)*kmrgd(3)-k3*kmrgd(1)*kmrgd(2)*kmrgd(4)+k8;
k1*kmrgd(1)*kmrgd(2)*kmrgd(3)-2*k2*kmrgd(3)^2+2*k3*kmrgd(1)*kmrgd(2)*kmrgd(4)-k4*kmrgd(3)+k6;
-k3*kmrgd(1)*kmrgd(2)*kmrgd(4)+2*k2*kmrgd(3)^2-k5*kmrgd(4);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
y0=[0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
handles = feval(SteinmetzLarter);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
jac=[[-k1*kmrgd(2)*kmrgd(3)-k3*kmrgd(2)*kmrgd(4)-km7,-k1*kmrgd(1)*kmrgd(3)-k3*kmrgd(1)*kmrgd(4),-k1*kmrgd(1)*kmrgd(2),-k3*kmrgd(1)*kmrgd(2)];[-k1*kmrgd(2)*kmrgd(3)-k3*kmrgd(2)*kmrgd(4),-k1*kmrgd(1)*kmrgd(3)-k3*kmrgd(1)*kmrgd(4),-k1*kmrgd(1)*kmrgd(2),-k3*kmrgd(1)*kmrgd(2)];[k1*kmrgd(2)*kmrgd(3)+2*k3*kmrgd(2)*kmrgd(4),k1*kmrgd(1)*kmrgd(3)+2*k3*kmrgd(1)*kmrgd(4),k1*kmrgd(1)*kmrgd(2)-4*k2*kmrgd(3)-k4,2*k3*kmrgd(1)*kmrgd(2)];[-k3*kmrgd(2)*kmrgd(4),-k3*kmrgd(1)*kmrgd(4),4*k2*kmrgd(3),-k3*kmrgd(1)*kmrgd(2)-k5]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
jacp=[[-kmrgd(1)*kmrgd(2)*kmrgd(3),0,-kmrgd(1)*kmrgd(2)*kmrgd(4),0,0,0,1,0,-kmrgd(1)];[-kmrgd(1)*kmrgd(2)*kmrgd(3),0,-kmrgd(1)*kmrgd(2)*kmrgd(4),0,0,0,0,1,0];[kmrgd(1)*kmrgd(2)*kmrgd(3),-2*kmrgd(3)^2,2*kmrgd(1)*kmrgd(2)*kmrgd(4),-kmrgd(3),0,1,0,0,0];[0,2*kmrgd(3)^2,-kmrgd(1)*kmrgd(2)*kmrgd(4),0,-kmrgd(4),0,0,0,0]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
hess1=[[0,-k1*kmrgd(3)-k3*kmrgd(4),-k1*kmrgd(2),-k3*kmrgd(2)];[0,-k1*kmrgd(3)-k3*kmrgd(4),-k1*kmrgd(2),-k3*kmrgd(2)];[0,k1*kmrgd(3)+2*k3*kmrgd(4),k1*kmrgd(2),2*k3*kmrgd(2)];[0,-k3*kmrgd(4),0,-k3*kmrgd(2)]];
hess2=[[-k1*kmrgd(3)-k3*kmrgd(4),0,-k1*kmrgd(1),-k3*kmrgd(1)];[-k1*kmrgd(3)-k3*kmrgd(4),0,-k1*kmrgd(1),-k3*kmrgd(1)];[k1*kmrgd(3)+2*k3*kmrgd(4),0,k1*kmrgd(1),2*k3*kmrgd(1)];[-k3*kmrgd(4),0,0,-k3*kmrgd(1)]];
hess3=[[-k1*kmrgd(2),-k1*kmrgd(1),0,0];[-k1*kmrgd(2),-k1*kmrgd(1),0,0];[k1*kmrgd(2),k1*kmrgd(1),-4*k2,0];[0,0,4*k2,0]];
hess4=[[-k3*kmrgd(2),-k3*kmrgd(1),0,0];[-k3*kmrgd(2),-k3*kmrgd(1),0,0];[2*k3*kmrgd(2),2*k3*kmrgd(1),0,0];[-k3*kmrgd(2),-k3*kmrgd(1),0,0]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
hessp1=[[-kmrgd(2)*kmrgd(3),-kmrgd(1)*kmrgd(3),-kmrgd(1)*kmrgd(2),0];[-kmrgd(2)*kmrgd(3),-kmrgd(1)*kmrgd(3),-kmrgd(1)*kmrgd(2),0];[kmrgd(2)*kmrgd(3),kmrgd(1)*kmrgd(3),kmrgd(1)*kmrgd(2),0];[0,0,0,0]];
hessp2=[[0,0,0,0];[0,0,0,0];[0,0,-4*kmrgd(3),0];[0,0,4*kmrgd(3),0]];
hessp3=[[-kmrgd(2)*kmrgd(4),-kmrgd(1)*kmrgd(4),0,-kmrgd(1)*kmrgd(2)];[-kmrgd(2)*kmrgd(4),-kmrgd(1)*kmrgd(4),0,-kmrgd(1)*kmrgd(2)];[2*kmrgd(2)*kmrgd(4),2*kmrgd(1)*kmrgd(4),0,2*kmrgd(1)*kmrgd(2)];[-kmrgd(2)*kmrgd(4),-kmrgd(1)*kmrgd(4),0,-kmrgd(1)*kmrgd(2)]];
hessp4=[[0,0,0,0];[0,0,0,0];[0,0,-1,0];[0,0,0,0]];
hessp5=[[0,0,0,0];[0,0,0,0];[0,0,0,0];[0,0,0,-1]];
hessp6=[[0,0,0,0];[0,0,0,0];[0,0,0,0];[0,0,0,0]];
hessp7=[[0,0,0,0];[0,0,0,0];[0,0,0,0];[0,0,0,0]];
hessp8=[[0,0,0,0];[0,0,0,0];[0,0,0,0];[0,0,0,0]];
hessp9=[[-1,0,0,0];[0,0,0,0];[0,0,0,0];[0,0,0,0]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
function userfun1=test1(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,k8,km7)
	userfun1=kmrgd(3)-0.018;
