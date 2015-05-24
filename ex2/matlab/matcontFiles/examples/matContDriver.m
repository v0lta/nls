clear all

%addpath(genpath('/home/r0483530/nonLinearSys/matlab/matcont6p1'))
addpath(genpath('/home/moritz/uni/nls/matlab/cl_matcont2.5.1'))
%addpath(genpath('/home/moritz/uni/nls/matlab/matcont5p3'))
%addpath(genpath('/home/moritz/uni/nls/matlab/matcont5p3'))
%addpath(genpath('/home/moritz/uni/nls/matlab/matcont2.4.2'))

global cds sys

%%%%% Set continuation pause environment variables %%%%%
%%
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

%%%%% Set system %%%%%
syshandle=@pitchfork;  %Specify system file


%%
%%%%% ODE Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Define a few intermediate functions for ODE integration %%%%%
SubFunHandles=feval(syshandle);  %Get function handles from system file
RHShandle=SubFunHandles{2};      %Get function handle for ODE

r=-1; %Set ODE parameter
xinit=1; %Set ODE initial condition

%%%%%% Define an anynomous function to pass to integrator %%%%%
RHS_no_param=@(t,x)RHShandle(t,x,r); 

%%%%% Set ODE integrator options %%%%%
options=odeset;
options=odeset(options,'RelTol',1e-5);
options=odeset(options,'maxstep',1e-1);

%%%%%% Integrate until a steady state is found. %%%%%
[tout xout]=ode45(RHS_no_param,[0,100],xinit,options);

figure()
plot(tout,xout,'-k','linewidth',2)
axis([0 100 -1 1])

%%
%%%%% Continuation from equilibrium %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Set initial condition as the endpoint of integration.  Use
%%%%%% to bootstrap the continuation.
xinit=xout(length(xout));


pvec=r; % Initialize parameter vector

ap=1; % Denote 'active' parameter for continuation

%%%%% Initialize continuation %%%%%
[x0,v0]=init_EP_EP(syshandle, xinit, pvec, ap); %Initialize equilibrium

%%%%% Initialize Matcont options %%%%%
opt=contset;
opt=contset(opt,'MaxNumPoints',200); %Set numeber of continuation steps
opt=contset(opt,'MaxStepsize',.01);  %Set max step size
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues 
opt = contset(opt,'InitStepsize',0.01); %Set Initial stepsize

%%%%% Continuation %%%%%
[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt); %Equilibrium continuation

% [x1,v1,s1,h1,f1]=cont(x1,v1,s1,h1,f1,cds);

figure()
cpl(x1,v1,s1,[2 1]);

%%
%%%%% Branch swiching anc continuation %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xBP=x1(1,s1(2).index);       %Extract branch point
pvec(ap)=x1(2,s1(2).index);  %Extract branch point 

%%%%%% Initialize branch point information %%%%%
[x0,vO]=init_BP_EP(syshandle, xBP, pvec, s1(2), 0.01);  


[x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt); %Switch branches and continue.

%figure()
%cpl(x2,v2,s2,[2 1]);

%%
%%%%% Backward continuation from near BP %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pick a point on the upper branch and use it to continue backward.
[x0,v0]=init_EP_EP(syshandle, x2(1,10), x2(2,10), ap);

% [x3,v3,s3,h3,f3]=cont(@equilibrium,x0,[],opt);

opt=contset(opt,'Backward',1);  %Compute backward
[x3,v3,s3,h3,f3]=cont(@equilibrium,x0,[],opt);

%figure()
%cpl(x3,v3,s3,[2 1]);

%%
%%%%% Plotting script.  x=continuation info  f=eigenvalues %%%%%%%%%
%%%%% Extract eigenvalues and match with curves for stability %%%%%%

figure()

xeqcurve=x1;

%%%%% This is the last eigenvalue.  
%%%%% These are ordered smallest to largest (real part).
%%%%% So the last one determines stability.
minevaleq=f1(1,:); 

L=length(xeqcurve(1,:));

curveind=1;
lengthind=0;
maxlengthind=0;
evalstart=floor(heaviside(minevaleq(1)));
datamateq=zeros(4,L);

for i=1:L
    evalind=floor(heaviside(minevaleq(i)));
    if evalstart~=evalind
        curveind=curveind+1;
        i;
        evalstart=evalind;
        maxlengthind=max(lengthind,maxlengthind);   
        lengthind=0;
    end
    datamateq(1,i)=xeqcurve(2,i); % This is the parameter that is varied.
    datamateq(2,i)=xeqcurve(1,i); % This is the dependent axis of the bifurcation plot.  The one you wish to plot
    datamateq(3,i)=evalind;
    datamateq(4,i)=curveind;  
    
    lengthind=lengthind+1;
end

maxlengthind=max(maxlengthind,lengthind);

curveindeq=curveind;

for i=1:curveindeq
    index=find(datamateq(4,:)==i);
    eval(['curve' num2str(i) 'eq' '=datamateq(1:3,index);']);
end


for i=1:curveindeq
    stability=eval(['curve' num2str(i) 'eq(3,1)']);
    if stability==0
        plotsty='-';
    else
        plotsty=':';
    end
    
    plotcolor='k';

    plotstr=strcat(plotcolor,plotsty);
    
    plot(eval(['curve' num2str(i) 'eq(1,:)']),eval(['curve' num2str(i) 'eq(2,:)']),plotstr','Linewidth',4)
    hold on
end


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
xeqcurve=x2;
minevaleq=f2(1,:); %This is the last eigenvalue.  That is the one that determines stability

L=length(xeqcurve(1,:));

curveind=1;
lengthind=0;
maxlengthind=0;
evalstart=floor(heaviside(minevaleq(1)));
datamateq=zeros(4,L);

for i=1:L
    evalind=floor(heaviside(minevaleq(i)));
    if evalstart~=evalind
        curveind=curveind+1;
        i;
        evalstart=evalind;
        maxlengthind=max(lengthind,maxlengthind);   
        lengthind=0;
    end
    datamateq(1,i)=xeqcurve(2,i); % This is the parameter that is varied.
    datamateq(2,i)=xeqcurve(1,i); % This is the dependent axis of the bifurcation plot.  The one you wish to plot
    datamateq(3,i)=evalind;
    datamateq(4,i)=curveind;  
    
    lengthind=lengthind+1;
end

maxlengthind=max(maxlengthind,lengthind);

curveindeq=curveind;

for i=1:curveindeq
    index=find(datamateq(4,:)==i);
    eval(['curve' num2str(i) 'eq' '=datamateq(1:3,index);']);
end

for i=1:curveindeq
    stability=eval(['curve' num2str(i) 'eq(3,1)']);
    if stability==0
        plotsty='-';
    else
        plotsty=':';
    end
    
    plotcolor='k';

    plotstr=strcat(plotcolor,plotsty);
    
    plot(eval(['curve' num2str(i) 'eq(1,:)']),eval(['curve' num2str(i) 'eq(2,:)']),plotstr','Linewidth',4)
    hold on
end



%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
xeqcurve=x3;
minevaleq=f3(1,:); %This is the last eigenvalue.  That is the one that determines stability

L=length(xeqcurve(1,:));

curveind=1;
lengthind=0;
maxlengthind=0;
evalstart=floor(heaviside(minevaleq(1)));
datamateq=zeros(4,L);

for i=1:L
    evalind=floor(heaviside(minevaleq(i)));
    if evalstart~=evalind
        curveind=curveind+1;
        i;
        evalstart=evalind;
        maxlengthind=max(lengthind,maxlengthind);   
        lengthind=0;
    end
    datamateq(1,i)=xeqcurve(2,i); % This is the parameter that is varied.
    datamateq(2,i)=xeqcurve(1,i); % This is the dependent axis of the bifurcation plot.  The one you wish to plot
    datamateq(3,i)=evalind;
    datamateq(4,i)=curveind;  
    
    lengthind=lengthind+1;
end

maxlengthind=max(maxlengthind,lengthind);

curveindeq=curveind;

for i=1:curveindeq
    index=find(datamateq(4,:)==i);
    eval(['curve' num2str(i) 'eq' '=datamateq(1:3,index);']);
end

for i=1:curveindeq
    stability=eval(['curve' num2str(i) 'eq(3,1)']);
    if stability==0
        plotsty='-';
    else
        plotsty=':';
    end
    
    plotcolor='k';

    plotstr=strcat(plotcolor,plotsty);
    
    plot(eval(['curve' num2str(i) 'eq(1,:)']),eval(['curve' num2str(i) 'eq(2,:)']),plotstr','Linewidth',4)
    hold on
end

%%%%% Final adjustments to make things look nicer %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([-1,1,-1.5,1.5])
scalefactor=3;
fontsizevar=11;
fontsizevar=fontsizevar*scalefactor;

xlabel('r','fontsize',fontsizevar)
ylabel('x','fontsize',fontsizevar)
title('Pitchfork','fontsize',fontsizevar)
set(gcf,'units','inches','pos',[0 0 3.1*scalefactor 2.61*scalefactor])
set(gca,'fontsize',fontsizevar)




