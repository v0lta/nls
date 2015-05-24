% This files runs mathcont!

%addpath('/home/r0483530/nonLinearSys/matlab/matcont6p1/Equilibrium')
addpath(genpath('/home/r0483530/nonLinearSys/matlab/matcont6p1'))

%% My Stuff
y0 = [0.5];
p0 = [0];
[x0,v0] = init_EP_EP(@bridge, y0, p0, 1);


