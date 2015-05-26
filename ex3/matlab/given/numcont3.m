% Nonlinear Systems - 2013
% Example Script - Imperfect bifurcations - Numerical Continuation
%   Exercise: 2.3

% Use: works only for the strategy used in numcont2

% Version: March, 2013
% nico.scheerlinck@cs.kuleuven.be

clear all
close all
global cds

% ---------------------------------------------------
  opt=contset;
  opt=contset(opt,'MaxNumPoints',50); 
  opt=contset(opt,'Singularities',1);
  opt=contset(opt,'Eigenvalues',1);
% ---------------------------------------------------------
  r  = 1.0;    % initial value active parameter r
% ---------------------------------------------------------
  h  = -0.1;    % fixed value ==> USER setting
% ---------------------------------------------------------
  p0 = [r;h]; 
  ap = 1;      % Note: index of continuation parameter = 1
% ---------------------------------------------------------

%  Corresponding initial values for u can be found using: 
%  >> roots([-0.5 0 r h])
  
figure(1);clf;hold on  
% ----------------------------------------------------------
% FIND fold (LP)
% Compute first solution branch containing fold (LP)
    u0 = 1.3613; % initial value for u    
    % Initialisation Step
      [x0,v0] = init_EP_EP(@myfuncont,u0,p0,ap);
      % Compute Backward
        opt=contset(opt,'Backward',1);
        [x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);
        cpl(x1,v1,s1,[2 1]);
        xlabel('r');ylabel('u')
        title(['h = ' num2str(h)])

% Define x and p0 for detected fold point
% For implementation details see manual p. 57-60
  xnew=x1(1,s1(2).index);
  p0(ap)=x1(end,s1(2).index);

% Initialise continuation LP-LP
  [x0,v0]=init_LP_LP(@myfuncont,xnew,p0,[1 2]); 

% Start a backward and a forward fold continuation  
% from the first LP detected on the previous equilibrium equation
  opt=contset(opt,'Backward',1);  % backward
  [x2,v2,s2,h2,f2]=cont(@limitpoint,x0,v0,opt); 
  opt=contset(opt,'Backward',0);  % forward
  [x3,v3,s3,h3,f3]=cont(@limitpoint,x0,v0,opt);

% PROCESS RESULTS  
% (u,r) - plane
figure(2);clf;hold on
cpl(x2,v2,s2,[1,2]);
cpl(x3,v3,s3,[1,2]);
axis([-2 2 -0.5 1])
xlabel('u');ylabel('r')
grid on

% (u,h) - plane
figure(3);clf;hold on
cpl(x2,v2,s2,[1,3]);
cpl(x3,v3,s3,[1,3]);
axis([-2 2 -1 1])
xlabel('u');ylabel('h')
grid on

% (r,h) - plane
figure(4);clf;hold on
cpl(x2,v2,s2,[2,3]);
cpl(x3,v3,s3,[2,3]);
axis([-0.5 1 -1 1])
xlabel('r');ylabel('h')
grid on

% (r,h,u) - space
figure(5);clf;hold on
cpl(x2,v2,s2,[2,3,1]);
cpl(x3,v3,s3,[2,3,1]);
axis([-0.5 1 -1 1 -2 2])
xlabel('r');ylabel('h');zlabel('u')
grid on

% (u,r,h) - space
figure(6);clf;hold on
cpl(x2,v2,s2,[1,2,3]);
cpl(x3,v3,s3,[1,2,3]);
axis([-2 2 -0.5 1 -1 1])
xlabel('u');ylabel('r');zlabel('h')
grid on