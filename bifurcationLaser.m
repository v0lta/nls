clear;

%The type of the bifuraction...

G = 1;
f  = 1;
k = 1;
p = 2;

n = 0:0.01:5;

%Determine the values of the original equation for different n
for i = 1:length(n)
dn(i) = G*n(i)*(p/(G*n(i)+f)) - k*n(i);
end

%figure(1)
%plot(n,dn)
%xlabel('n')
%ylabel('dn')


pc = k*f/G;
p = 0:0.1:2*pc;

%Compute the values of the second derivative.
%figure(2)
for i = 1:length(p)
   nTwo = (G*p(i)-f*k)/k*G;
%   plot(p(i),nTwo,'b*')
%   hold on
end
%plot(0,0,'r*')


%find the zeros of the orginal equation.
%for different coeficients and plot them.
for i = 1:length(p)
 funDn = @(n) G*n*(p(i)/(G*n+f)) - k*n;
 zero(i) = fzero(funDn,pc);
 for j = 1:length(n)
    dn(i,j) = funDn(n(j));
 end
end

subplot(1,2,1)
plot(p,real(zero),'*');
xlabel('p')
ylabel('zero')

subplot(1,2,2)
plot(n,dn)
hold on;
plot(zero,0,'r*')
xlabel('n')
ylabel('dot n')


%Its a transcritical bifurcation.
