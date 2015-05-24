% This files runs mathcont!

%addpath('/home/r0483530/nonLinearSys/matlab/matcont6p1/Equilibrium')
addpath(genpath('/home/r0483530/nonLinearSys/matlab/matcont6p1'))

%% My Stuff
%y0 = [0.5 0]';
%p0 = [0 0]';
%x0,v0] = init_EP_EP(@bridge, y0, p0, 1);

%% The pitch example.
%setting some options see manual p. 14f
opt = contset;
opt = contset(opt, 'MaxNumPoints', 30);
opt = contset(opt, 'Singularities',1);

% run first continuation
[x0, v0] = cont(@pitch, [0], [-1], [1]);