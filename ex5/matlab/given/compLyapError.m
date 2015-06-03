


%x = [0 1 2500]';
%x = [0 1 20]';
%x = [0.5 0.5 0.5]';
%x = [0.1 0.1 0.1]';
x = [1 1 1]';




st = 0.001;
kkmax = 1000;
lyap = lyapunov(@rhs_lorenz,st,kkmax,x);
