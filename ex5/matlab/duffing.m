%% compute the data.
tend = 20;
[T1,X1] = ode45(@duffFun,[0 tend],[2.0 0]);
[T2,X2] = ode45(@duffFun,[0 tend],[2.1 0]);
[T3,X3] = ode45(@duffFun,[0 tend],[2.2 0]);
[T4,X4] = ode45(@duffFun,[0 tend],[2.3 0]);

%arrange some nice colors.
figure(1)
colorMat = lines;
color1 = colorMat(1,:);
color2 = colorMat(2,:);
color3 = colorMat(3,:);
color4 = colorMat(4,:);

plot(X1(:,1),X1(:,2),'color',color1)
hold on
plot(X2(:,1),X2(:,2),'color',color2)
plot(X3(:,1),X3(:,2),'color',color3)
plot(X4(:,1),X4(:,2),'color',color4)

plot(X1(end,1),X1(end,2),'o','color',color1)
plot(X2(end,1),X2(end,2),'o','color',color2)
plot(X3(end,1),X3(end,2),'o','color',color3)
plot(X4(end,1),X4(end,2),'o','color',color4)
hold off;

figure(2)
plot(T1,X1(:,1),'color',color1)
plot(T1,X1(:,2),'color',color1)
hold on
plot(T2,X2(:,1),'color',color2)
plot(T2,X2(:,2),'color',color2)
plot(T3,X3(:,1),'color',color3)
plot(T3,X3(:,2),'color',color3)
plot(T4,X4(:,1),'color',color4)
plot(T4,X4(:,2),'color',color4)


% %% make a movie
% writerObj = VideoWriter('duffing.avi');
% open(writerObj);
% 
% axis tight
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
% 
% for k = 1:(length(X4))
%    plot(X1(k,1),X1(k,2),'o','color',color1)
%    hold on;
%    plot(X2(k,1),X2(k,2),'o','color',color2)
%    plot(X3(k,1),X3(k,2),'o','color',color3)
%    plot(X4(k,1),X4(k,2),'o','color',color4)
%    frame = getframe;
%    writeVideo(writerObj,frame);
% end
% hold off;
% 
% close(writerObj);