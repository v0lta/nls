

G = 1;
k = 1;
f = 1;

p = -1:0.5:3;

size = 25;




for i=1:1:length(p)
    
       
    if (p(i) < (f*k/G))
        x1 = [ 0; p(i)/f];
        x2 = [ -f/G + p(i)/k ; k/G];
        plot(real(x1(1)),real(x1(2)),'b.','MarkerSize',size);
        hold on;
        plot(real(x2(1)),real(x2(2)),'ro');
    end 
    
    
    if (p(i) == (f*k/G))
        x1 = [ 0; p(i)/f];
        x2 = [ -f/G + p(i)/k ; k/G];
        plot(real(x1(1)),real(x1(2)),'k.','MarkerSize',size);
        hold on;
        plot(real(x2(1)),real(x2(2)),'k.','MarkerSize',size);
    end 
    
    if (p(i) > (f*k/G))
        x1 = [ 0; p(i)/f];
        x2 = [ -f/G + p(i)/k ; k/G];
        plot(real(x1(1)),real(x1(2)),'bo');
        hold on;
        plot(real(x2(1)),real(x2(2)),'r.','MarkerSize',size);
    end 
    
end
grid on;


