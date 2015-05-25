% This files runs mathcont!

addpath(genpath('../../../../matlab/cl_matcont2.5.1/'))

clear all
global cds

opt=contset;
opt=contset(opt,'MaxNumPoints',30);
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MinStepsize',0.1);
opt=contset(opt,'MaxStepsize',10.0);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);

%% My Stuff
y0 = [1 0]';
V0 = 42.5985/2;
[x0,v0] = init_EP_EP(@bridge, y0, V0, 1);
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);

cpl(x1,v1,s1,[1 2]);

%%  -----------------------------------------------
xh = x1(1:2,s1(2).index);
Vh = x1(3,s1(2).index);
i  = s1(2).index;
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

  [x0,v0]=init_H_LC(@bridge,xh,Vh,1,1e-6,ntst,ncol);
  [x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);

  
figure(2);clf;
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


