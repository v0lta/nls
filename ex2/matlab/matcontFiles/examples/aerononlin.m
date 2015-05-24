% Nonlinear Systems - 2013
% Example Script - Aero-elastic galloping - Numerical Continuation 
%   Exercise: 3.6
%   System:   see funaero.m in Systems folder
%             Note: the state variables y and ydot have been adjusted
%                   relative to Vc, i.e., the critical windspeed. 
%
%!! Goal:     Bifurcation diagram of limit cycles
%
%*****************************************************
%!! Use:      Example script -- NOT working properly 
%             ==> please refine and extend the script
%*****************************************************
%!! More details about MATCONT routines: 
%   see CL MATCONT manual
%       Topic: Continuation of limit cycles: p. 46-56

% Version: March, 2013
% nico.scheerlinck@cs.kuleuven.be

clear all
global cds

opt=contset;
opt=contset(opt,'MaxNumPoints',30);
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MinStepsize',0.1);
opt=contset(opt,'MaxStepsize',10.0);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);

disp('Computations are running ...')
disp('Please wait ...')
time = cputime;
p0  = 0.2;  % initial parameter value
[x0,v0]          = init_EP_EP(@funaero,[0.5;0],p0,1);
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);

xh = x1(1:2,s1(2).index);
Vh = x1(3,s1(2).index);
i  = s1(2).index;

% -----------------------------------------------
% IMPORTANT
  opt=contset(opt,'MinStepsize',0.001);
  opt=contset(opt,'Multipliers',1);
  %opt=contset(opt,'MaxNewtonIters',5);
  opt=contset(opt,'Adapt',1);
  opt=contset(opt,'MaxNumPoints',120);

  ntst = 8;  % number of mesh-intervals
             % this is an initial guess !!!!!!
             % increase this argument !
  ncol = 4;  % number of collocation points
% -----------------------------------------------
  disp('Computations are still running ...')
  disp('Please wait ...')

  [x0,v0]=init_H_LC(@funaero,xh,Vh,1,1e-6,ntst,ncol);
  [x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);

figure(1);clf;
[m points] = size(x2); % # of continuation steps
dim = 2; % dimension of ODE
% max x2 along periodic solutions
xx = x2(1:end-2,:);
xx = reshape(xx, [dim (m-2)/dim points]);
%xa = max(squeeze(xx(1,:,:)), [], 1)/42.6;
xa = max(squeeze(xx(1,:,:)), [], 1);
%mu = x2(end,:)/42.6;
mu = x2(end,:);

cpl([mu; xa; zeros(size(xa))],v2,s2,[1 2]);
axis([0 2.5 -.5 4])

figure(2);clf;
e = size(x2,1)-1;
plotcycle(x2([1:end-2 end],1:end),v2,s2,[e 1 2])
view(3)
xlabel('V/Vc')
ylabel('y/Vc')
zlabel('ydot/Vc')
grid on

figure(3);clf;
e = size(x2,1)-1;
A1 = 0.04695;
Vc = 2/A1;
new = [x2([1:end-2],1:end)*Vc; x2(end,1:end)];
plotcycle(new,v2,s2,[e 1 2])
view(3)
xlabel('V/Vc')
ylabel('y')
zlabel('ydot')
grid on

timesec = cputime-time;
timemin=timesec/60
disp('Computations are DONE !')