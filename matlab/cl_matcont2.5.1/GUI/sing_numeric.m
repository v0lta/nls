function e=sing_numeric(ignore)
%de functie bepaalt de index in de h
%
len=size(ignore,2);
for i=1:8
    e(1,i)=2+i;
end
ignore=sort(ignore);
for i=1:len
    switch ignore(i)
    case 1
        e(1,1)=0;
        e(1,2)=e(1,2)-1;
        e(1,3)=e(1,3)-1;
        e(1,4)=e(1,4)-1;
        e(1,5)=e(1,5)-1;
        e(1,6)=e(1,6)-1;
        e(1,7)=e(1,7)-1;
        e(1,8)=e(1,8)-1;

    case 2
        e(1,2)=0;
        e(1,3)=e(1,3)-1;
        e(1,4)=e(1,4)-1;
        e(1,5)=e(1,5)-1;
        e(1,6)=e(1,6)-1;
        e(1,7)=e(1,7)-1;
        e(1,8)=e(1,8)-1;
    case 3
        e(1,3)=0;
        e(1,4)=e(1,4)-1;
        e(1,5)=e(1,5)-1;
        e(1,6)=e(1,6)-1;
        e(1,7)=e(1,7)-1;
        e(1,8)=e(1,8)-1;
    case 4
        e(1,4)=0;
        e(1,5)=e(1,5)-1;
        e(1,6)=e(1,6)-1;
        e(1,7)=e(1,7)-1;
        e(1,8)=e(1,8)-1;
    case 5
        e(1,5)=0;
        e(1,6)=e(1,6)-1;
        e(1,7)=e(1,7)-1;
        e(1,8)=e(1,8)-1;
    case 6
        e(1,6)=0;
        e(1,7)=e(1,7)-1;
        e(1,8)=e(1,8)-1;
    case 7
        e(1,7)=0;
        e(1,8)=e(1,8)-1;
    case 8
        e(1,8)=0;
    end
end
