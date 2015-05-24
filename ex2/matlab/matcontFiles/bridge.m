function [ out ] = bridge
%bridge starts the mathcont computations for problem two using the mathcont
%reference from schilda :-).

out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
end

%-------------------------------------------------------------------------
function dy = fun_eval(t,y,V) 

m = 1;
rho = 1;
r = 1;
k = 100;
a = 1;

A1 = 0.04695;
A3 = 8.932*10^(-4);
A5 = 1.015*10^(-5);
A7 = 2.955*10^(-8);

dy(1) = y(2);
dy(2) = -k/m * y(1) - r/m * y(2) + (0.5/m)*(rho*V^2*a) ...
            * (A1 * (y(2)/V)^1 - A3 * (y(2)/V)^3 + A5 * (y(2)/V)^5 ...
            - A7 * (y(2)/V)^7);
end

%-------------------------------------------------------------------------
function [tspan,y0,options] = init
y0 = [0,0];
options = odeset;
handles = feval(bridge);
tspan = [0 10];
end
