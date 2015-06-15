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
opt=contset(opt,'MinStepsize',0.01);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);

disp('Computations are running ...')
disp('Please wait ...')
time = cputime;

c = 0.1;
%d = 0;
d = +0.01;
%d = -0.01;
plot = [2 3];

c0  = [c; d];  % initial parameter value

%x4 [0.4; 0]
if (d == 0)
    opt=contset(opt,'MaxNumPoints',25);
else
    opt=contset(opt,'MaxNumPoints',50);
end
[x0,v0]          = init_EP_EP(@funPred,[0.4;0],c0,1);
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);
[x0,v0]          = init_EP_EP(@funPred,[0.4;0],[1; d],1);
[x1b,v1b,s1b,h1b,f1b] = cont(@equilibrium,x0,v0,opt);

%x3 [1; 0]

opt=contset(opt,'MaxNumPoints',30);
[x0,v0]          = init_EP_EP(@funPred,[1;0],c0,1);
[x3,v3,s3,h3,f3] = cont(@equilibrium,x0,v0,opt);
[x0,v0]          = init_EP_EP(@funPred,[1.0;0],[1.25; d],1);
opt=contset(opt,'Backward',1);
[x3b,v3b,s3b,h3b,f3b] = cont(@equilibrium,x0,v0,opt);
opt=contset(opt,'Backward',0);
[x3c,v3c,s3c,h3c,f3c] = cont(@equilibrium,x0,v0,opt);
%x2 
%from bottom.
if (d == 0) 
    opt=contset(opt,'MaxNumPoints',50);
else
    opt=contset(opt,'MaxNumPoints',50);
end
[x0,v0]          = init_EP_EP(@funPred,[0.1;-0.9],c0,1);
[x2,v2,s2,h2,f2] = cont(@equilibrium,x0,v0,opt);

%split approach from the top.
u0 = [0.7; 0.3]; % initial value for u    
%   %Initialisation Step
[x0,c0] = init_EP_EP(@funPred,u0,c0,1);
%   % Compute Forward
opt=contset(opt,'Backward',0);
[x2a,v2a,s2a,h2a,f2a] = cont(@equilibrium,x0,c0,opt);
%   % Compute Backward
opt=contset(opt,'Backward',1);
[x2b,v2b,s2b,h2b,f2b] = cont(@equilibrium,x0,c0,opt);
   

figure(4)
if (d == 0)
    cpl(x2,v2,s2,plot);
else
    cpl(x2,v2,s2,plot);
    cpl(x2a,v2a,s2a,plot);
    cpl(x2b,v2b,s2b,plot);
end
hold on;
%cpl(x1,v1,s1,plot);
%cpl(x3,v3,s3,plot);
cpl(x3b,v3b,s3b,plot)
cpl(x3c,v3c,s3c,plot)
%cpl(x1b,v1b,s1b,plot)
title(['d = ' num2str(d)])
hold off;

%select values for limit cycle continuation.
if (d == 0)
    xh = x2(1:2,s2(3).index);
    ch = [x2(3,s2(3).index); d];
    i  = s2(2).index;
elseif (d == +0.01)
    xh = x2a(1:2,s2a(2).index);
    ch = [x2a(3,s2a(2).index); d];
else
    xh = x2b(1:2,s2b(2).index);
    ch = [x2b(3,s2b(2).index); d];
end

%% Part2 -----------------------------------------------
% IMPORTANT
  opt=contset(opt,'MinStepsize',0.01);
  opt=contset(opt,'Multipliers',1);
  %opt=contset(opt,'MaxNewtonIters',5);
  opt=contset(opt,'Adapt',1);
  opt=contset(opt,'MaxNumPoints',120);


  ntst = 15;  % number of mesh-intervals
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