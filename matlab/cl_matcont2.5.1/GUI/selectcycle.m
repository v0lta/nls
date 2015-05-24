function [x,v,s,h,f] = selectcycle(varargin)
% "Grab" periodic orbit after time integration, to start LC continuation.
% Syntax:
%   [x,v,s,h,f] = selectcycle(t, x, tolerance, ntst, param)
% x should contain data from time integration, where it STARTED ON the
% orbit, and does NOT circle the orbit TWICE OR MORE.

global lds
t = varargin{1};
x = varargin{2}';
tolerance = varargin{3};
ntst = varargin{4};
param = varargin{5};
xstart=round(size(x,2)/3);
amin=sum(abs(x(:,xstart:end)-x(:,1)*ones(1,size(x(:,xstart:end),2))));
ep=tolerance;
amin(find(amin(:)<ep))=inf;
[pp,qq]=min(amin-(1e-4));
if (pp==Inf)||pp<0
   ind=[]; 
else
    ind = qq(1)+xstart-1;
end
if isempty(ind)
    warndlg('No cycle can be found!')
    return; 
end
x=x(:,1:ind);
t=t(1:ind)';
tn=(t-t(1))/(t(end)-t(1));
[a,x,tn] = newmeshcycle(x,tn,size(x,2)-1,1,ntst,4);
x = interp(tn,1,x,a,4);
s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = 4;
s(1).data.parametervalues = param;
lds.P0 = param;
lds.ActiveParams = 1;
lds.ntst = ntst;
lds.ncol = 4;
x = reshape(x,size(x,2)*size(x,1),1);
x(end+1) = t(end)-t(1);
s(1).data.T = x(end);
lds.T = x(end);
x(end+1) = param(1);
if size(s) == 2
    s(2)=[];
end
v=[];h=[];f=[];cds=[];ctype='LC ';point='LC ';num=1;
s.index = 1;
s(num).label = ctype;
