

% ---------------------------------------------------
  opt=contset;
  opt=contset(opt,'MaxNumPoints',500); 
  opt=contset(opt,'Singularities',1);
  opt=contset(opt,'Eigenvalues',1);
% ---------------------------------------------------------
  d  = 0;    % initial value active parameter r
% ---------------------------------------------------------
  c  = 0.7;    % fixed value ==> USER setting
% ---------------------------------------------------------
  p0 = [c;d]; 
  ap = 1;      % Note: index of continuation parameter = 1
% ---------------------------------------------------------


figure(1);hold on  
 % Compute first solution branch containing fold (LP)
 u0 = [0.7; 0.3]; % initial value for u    
 % Initialisation Step
 [x0,v0] = init_EP_EP(@myfuncont,u0,p0,ap);
 [x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);