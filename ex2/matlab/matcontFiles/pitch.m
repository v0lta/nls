function [ output_args ] = pitch
%The pitchfork example form the reference.


out{1} = @init;
out{2} = @funEval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
end

%-------------------------------------------------------------------------
function dydt = funEval(t,x,mu) 
dydt(1,1) = mu*x - x^3;
end

%-------------------------------------------------------------------------
function [tspan,y0,options] = init
y0 = [0,0];
options = odeset;
handles = feval(pitch);
tspan = [0 10];
end


