

%The type of the bifuraction...

G = 1;
f  = 1;
k = 2;
p = 2;

n = -5:0.1:5;

for i = 1:length(n)
dn(i) = G*n(i)*(p/(G*n(i)+f)) - k*n(i);
end

figure(1)
plot(n,dn)

pc = k*f/G
p = 0:0.1:2*pc;

figure(2)
for i = 1:length(p)
   nTwo = (G*p(i)-f*k)/k*G
   plot(nTwo,0,'b*')
   hold on
end
plot(0,0,'r*')

%Its a transcritical bifurcation.
