% This files runs mathcont!

%addpath('/home/r0483530/nonLinearSys/matlab/matcont6p1/Equilibrium')
addpath(genpath('/home/r0483530/nonLinearSys/matlab/matcont6p1'))


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
y0 = [0.5 0]';
V0 = 42.5985*0.2;
[x0,v0] = init_EP_EP(@bridge, y0, V0, 1);
[x1,v1,s1,h1,f1] = cont(@equilibrium,x0,v0,opt);


