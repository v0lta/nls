% Nonlinear Systems - 2013
% Example Script - Imperfect bifurcations - Numerical Continuation
%   Exercise: 2.2

%   System:   Simplified equilibrium equation: 
%             u' = -0.5*u^3 + r*u + h = 0
%             see myfuncont.m in Systems folder
%!! Case:     r is used as a parameter, h is fixed 
%
%!! Goal:     Draw solution branches using r as parameter
%             for different values of h. 
%!! How?      This script has been written based on a first analysis
%             of the problem using MAPLE (see file on TOLEDO)
%!! Use:      This script can be used to draw solution branches
%             for h = 0, h > 0, h < 0
%
%!! More details about MATCONT routines: 
%   see CL MATCONT manual
%       Topic: Equilibrium continuation: p. 34-40

% Use: the code below has been tested for the following 3 case studies
%      1) h = -0.1 (fixed parameter) & r is the active parameter
%      2) h =  0   (fixed parameter) & r is the active parameter
%      3) h = +0.1 (fixed parameter) & r is the active parameter

%  Corresponding initial values for u can be found using: 
%  >> roots([-0.5 0 r h])

% Version: March, 2013
% nico.scheerlinck@cs.kuleuven.be

clear all
global cds

% ---------------------------------------------------------
% This is a first trial
% Adjust   'MaxNumPoints' if necessary 

  opt=contset;
  opt=contset(opt,'MaxNumPoints',400); 
  opt=contset(opt,'Singularities',1);
  opt=contset(opt,'Eigenvalues',1);
  % my stuff.
  opt=contset(opt,'MinStepsize',0.1);
  opt=contset(opt,'MaxStepsize',0.01);
  %watch the step stepsize ratio, parameter ratio.
  
  
  opt=contset(opt,'MaxNewtonIters',4);

  c = 0.7;
  d = 0.01;
  p0 = [c; d];
  ap = 1;      % Note: index of continuation parameter = 1
  plotVec = [1 2];
% ---------------------------------------------------------
% figure(1);clf;hold on  
% % Compute first solution branch
   u0 = [0.7; 0.3]; % initial value for u    
%   %Initialisation Step
     [x0,c0] = init_EP_EP(@myfuncont,u0,p0,ap);
%   % Compute Forward
     opt=contset(opt,'Backward',0);
     [x1,v1,s1,h1,f1] = cont(@equilibrium,x0,c0,opt);
     cpl(x1,v1,s1,plotVec);
   u0 = [0.7; 0.3]; % initial value for u    
%   %Initialisation Step
     [x0,c0] = init_EP_EP(@myfuncont,u0,p0,ap);
%   % Compute Forward
    opt=contset(opt,'Backward',1);
     [x1,v1,s1,h1,f1] = cont(@equilibrium,x0,c0,opt);
     cpl(x1,v1,s1,plotVec); 
     
