% Nonlinear Systems - 2013
% Example Script - Predetor prey model - Numerical Continuation 
%   Exercise: 3.6
%   System:   see funPred.m in Systems folder
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


%% Part1 -----------------------------------------------

clear all
global cds

opt=contset;
opt=contset(opt,'MaxNumPoints',25);
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MinStepsize',0.1);
opt=contset(opt,'MaxStepsize',1.0);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);

disp('Computations are running ...')
disp('Please wait ...')
time = cputime;
c0  = 0.1;  % initial parameter value
[x0,v0]          = init_EP_EP(@funPred,[0.1;-0.9],c0,1);
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);

xh = x1(1:2,s1(3).index);
ch = x1(3,s1(3).index);
i  = s1(2).index;

figure(4)
cpl(x1,v1,s1,[1 2]);

%% Part2 -----------------------------------------------
% IMPORTANT
  opt=contset(opt,'MinStepsize',0.01);
  opt=contset(opt,'Multipliers',1);
  %opt=contset(opt,'MaxNewtonIters',5);
  opt=contset(opt,'Adapt',1);
  opt=contset(opt,'MaxNumPoints',120);


  ntst = 20;  % number of mesh-intervals
             % this is an initial guess !!!!!!
             % increase this argument !
  ncol = 4;  % number of collocation points
% -----------------------------------------------
  disp('Computations are still running ...')
  disp('Please wait ...')

  [x0,v0]=init_H_LC(@funPred,xh,ch,1,1e-6,ntst,ncol);
  [x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);

figure(2);clf;
e = size(x2,1)-1;
plotcycle(x2([1:end-2 end],1:end),v2,s2,[e 1 2])
view(3)
xlabel('c')
ylabel('x')
zlabel('y')
grid on


timesec = cputime-time;
timemin=timesec/60
disp('Computations are DONE !')