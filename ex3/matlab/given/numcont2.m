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
  opt=contset(opt,'MaxNumPoints',500); 
  opt=contset(opt,'Singularities',1);
  opt=contset(opt,'Eigenvalues',1);
  % my stuff.
  opt=contset(opt,'MinStepsize',0.001);
  opt=contset(opt,'MaxStepsize',0.01);
  %watch the step stepsize ratio, parameter ratio.
  
  
  opt=contset(opt,'MaxNewtonIters',4);
% ---------------------------------------------------------
  r  = 1;    % initial value active parameter r
% ---------------------------------------------------------
  h  = -0.0005;    % fixed value ==> USER setting
% ---------------------------------------------------------
  p0 = [r;h]; 
  ap = 1;      % Note: index of continuation parameter = 1
% ---------------------------------------------------------
figure(1);clf;hold on  
% Compute first solution branch  
  u0 = 1.3613; % initial value for u    
  % Initialisation Step
    [x0,v0] = init_EP_EP(@myfuncont,u0,p0,ap);
  % Compute Backward
    opt=contset(opt,'Backward',1);
    [x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);
    cpl(x1,v1,s1,[2 1]);
% Compute remaining solution branch          
  u0 = -1.4618; % initial value for u    
  % Initialisation Step
    [x0,v0] = init_EP_EP(@myfuncont,u0,p0,ap);
  % Compute Forward
    opt=contset(opt,'Backward',0);
    [x2,v2,s2,h2,f2] = cont(@equilibrium,x0,v0,opt);
    cpl(x2,v2,s2,[2 1]);
    grid on
    xlabel('r');ylabel('u')
    title(['h = ' num2str(h)])

% POST PROCESSING AREA
  x = [x1 x2]; % forward and backward results
  f = [f1 f2]; % eigenvalues jacobian
  u = x(1,:);  % fixed points
  r = x(2,:);  % active parameter
  
  % Stability analysis plot
  % Goal: Visualisation of stable and unstable branches
    figure(2);clf;hold on
    pos = find(sign(f) == 1);   
    nul = find(sign(f) == 0);   
    neg = find(sign(f) == -1);
    
    if ~isempty(pos)
        plot(r(pos),u(pos),'r.','linewidth',2)
    end
    if ~isempty(nul)
        plot(r(nul),u(nul),'k.','linewidth',2)
    end
    if ~isempty(neg)
        plot(r(neg),u(neg),'b.','linewidth',2)
    end
    xlabel('r');ylabel('u')
    title(['h = ' num2str(h)])
    grid on
