function out = Ermentrout
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
alphaM=0.32*(54+kmrgd(1))/(1-exp(-(kmrgd(1)+54)/4));
alphaN=0.032*(kmrgd(1)+52)/(1-exp(-(kmrgd(1)+52)/5));
alphaH=0.128*exp(-(50+kmrgd(1))/18);
betaM=0.28*(kmrgd(1)+27)/(exp((kmrgd(1)+27)/5)-1);
betaN=0.5*exp(-(57+kmrgd(1))/40);
 betaH=4/(1+exp(-(kmrgd(1)+27)/5));
tw=100/(3.3*exp((kmrgd(1)+35)/20)+exp(-(kmrgd(1)+35)/20));
winf=1/(1+exp(-(kmrgd(1)+35)/10));
mlinf=1/(1+exp(-(kmrgd(1)+25)/2.5));
pca=gca*mlinf*(kmrgd(1)-eca);
GTOT=gna*kmrgd(2)^3*kmrgd(4)*(kmrgd(1)-ena)+(gk*kmrgd(3)^4+gm*kmrgd(5)+gahp*kmrgd(6)/(1+kmrgd(6)))*(kmrgd(1)-ek)+gl*(kmrgd(1)-el)+pca;
dydt=[(I-GTOT)/C;
(1-kmrgd(2))*alphaM-kmrgd(2)*betaM;
(1-kmrgd(3))*alphaN-kmrgd(3)*betaN;
(1-kmrgd(4))*alphaH-kmrgd(4)*betaH;
(winf-kmrgd(5))/tw;
-0.002*pca-kmrgd(6)/80;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Ermentrout);
y0=[0,0,0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
jac=[[(-gna*kmrgd(2)^3*kmrgd(4)-gk*kmrgd(3)^4-gm*kmrgd(5)-gahp*kmrgd(6)/(1+kmrgd(6))-gl-2/5*gca/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)-gca/(1+exp(-2/5*kmrgd(1)-10)))/-C,-3*gna*kmrgd(2)^2*kmrgd(4)*(kmrgd(1)-ena)/-C,-4*gk*kmrgd(3)^3*(kmrgd(1)-ek)/-C,-gna*kmrgd(2)^3*(kmrgd(1)-ena)/-C,-gm*(kmrgd(1)-ek)/-C,-(gahp/(1+kmrgd(6))-gahp*kmrgd(6)/(1+kmrgd(6))^2)*(kmrgd(1)-ek)/-C];[8/25*(1-kmrgd(2))/(1-exp(-27/2-1/4*kmrgd(1)))-1/4*(1-kmrgd(2))*(432/25+8/25*kmrgd(1))/(1-exp(-27/2-1/4*kmrgd(1)))^2*exp(-27/2-1/4*kmrgd(1))-7/25*kmrgd(2)/(exp(1/5*kmrgd(1)+27/5)-1)+1/5*kmrgd(2)*(7/25*kmrgd(1)+189/25)/(exp(1/5*kmrgd(1)+27/5)-1)^2*exp(1/5*kmrgd(1)+27/5),-(432/25+8/25*kmrgd(1))/(1-exp(-27/2-1/4*kmrgd(1)))-(7/25*kmrgd(1)+189/25)/(exp(1/5*kmrgd(1)+27/5)-1),0,0,0,0];[4/125*(1-kmrgd(3))/(1-exp(-1/5*kmrgd(1)-52/5))-1/5*(1-kmrgd(3))*(4/125*kmrgd(1)+208/125)/(1-exp(-1/5*kmrgd(1)-52/5))^2*exp(-1/5*kmrgd(1)-52/5)+1/80*kmrgd(3)*exp(-57/40-1/40*kmrgd(1)),0,-(4/125*kmrgd(1)+208/125)/(1-exp(-1/5*kmrgd(1)-52/5))-1/2*exp(-57/40-1/40*kmrgd(1)),0,0,0];[-8/1125*(1-kmrgd(4))*exp(-25/9-1/18*kmrgd(1))-4/5*kmrgd(4)/(1+exp(-1/5*kmrgd(1)-27/5))^2*exp(-1/5*kmrgd(1)-27/5),0,0,-16/125*exp(-25/9-1/18*kmrgd(1))-4/(1+exp(-1/5*kmrgd(1)-27/5)),0,0];[1/1000/(1+exp(-1/10*kmrgd(1)-7/2))^2*exp(-1/10*kmrgd(1)-7/2)*(33/10*exp(1/20*kmrgd(1)+7/4)+exp(-1/20*kmrgd(1)-7/4))+1/100*(1/(1+exp(-1/10*kmrgd(1)-7/2))-kmrgd(5))*(33/200*exp(1/20*kmrgd(1)+7/4)-1/20*exp(-1/20*kmrgd(1)-7/4)),0,0,0,-33/1000*exp(1/20*kmrgd(1)+7/4)-1/100*exp(-1/20*kmrgd(1)-7/4),0];[-1/1250*gca/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)-1/500*gca/(1+exp(-2/5*kmrgd(1)-10)),0,0,0,0,-1/80]];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
jacp=[[-(I-gna*kmrgd(2)^3*kmrgd(4)*(kmrgd(1)-ena)-(gk*kmrgd(3)^4+gm*kmrgd(5)+gahp*kmrgd(6)/(1+kmrgd(6)))*(kmrgd(1)-ek)-gl*(kmrgd(1)-el)-gca/(1+exp(-2/5*kmrgd(1)-10))*(kmrgd(1)-eca))/-C^2,1/-C,-kmrgd(2)^3*kmrgd(4)*(kmrgd(1)-ena)/-C,gna*kmrgd(2)^3*kmrgd(4)/-C,-kmrgd(3)^4*(kmrgd(1)-ek)/-C,(gk*kmrgd(3)^4+gm*kmrgd(5)+gahp*kmrgd(6)/(1+kmrgd(6)))/-C,(-kmrgd(1)+el)/-C,gl/-C,-kmrgd(5)*(kmrgd(1)-ek)/-C,-kmrgd(6)/(1+kmrgd(6))*(kmrgd(1)-ek)/-C,-1/(1+exp(-2/5*kmrgd(1)-10))*(kmrgd(1)-eca)/-C,gca/(1+exp(-2/5*kmrgd(1)-10))/-C];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,0,0];[0,0,0,0,0,0,0,0,0,0,-1/500/(1+exp(-2/5*kmrgd(1)-10))*(kmrgd(1)-eca),1/500*gca/(1+exp(-2/5*kmrgd(1)-10))]];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
hess1=[[(-8/25*gca/(1+exp(-2/5*kmrgd(1)-10))^3*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)^2-4/5*gca/(1+exp(-2/5*kmrgd(1)-10))^2*exp(-2/5*kmrgd(1)-10)+4/25*gca/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10))/-C,-3*gna*kmrgd(2)^2*kmrgd(4)/-C,-4*gk*kmrgd(3)^3/-C,-gna*kmrgd(2)^3/-C,-gm/-C,-(gahp/(1+kmrgd(6))-gahp*kmrgd(6)/(1+kmrgd(6))^2)/-C];[-4/25*(1-kmrgd(2))/(1-exp(-27/2-1/4*kmrgd(1)))^2*exp(-27/2-1/4*kmrgd(1))+1/8*(1-kmrgd(2))*(432/25+8/25*kmrgd(1))/(1-exp(-27/2-1/4*kmrgd(1)))^3*exp(-27/2-1/4*kmrgd(1))^2+1/16*(1-kmrgd(2))*(432/25+8/25*kmrgd(1))/(1-exp(-27/2-1/4*kmrgd(1)))^2*exp(-27/2-1/4*kmrgd(1))+14/125*kmrgd(2)/(exp(1/5*kmrgd(1)+27/5)-1)^2*exp(1/5*kmrgd(1)+27/5)-2/25*kmrgd(2)*(7/25*kmrgd(1)+189/25)/(exp(1/5*kmrgd(1)+27/5)-1)^3*exp(1/5*kmrgd(1)+27/5)^2+1/25*kmrgd(2)*(7/25*kmrgd(1)+189/25)/(exp(1/5*kmrgd(1)+27/5)-1)^2*exp(1/5*kmrgd(1)+27/5),-8/25/(1-exp(-27/2-1/4*kmrgd(1)))+1/4*(432/25+8/25*kmrgd(1))/(1-exp(-27/2-1/4*kmrgd(1)))^2*exp(-27/2-1/4*kmrgd(1))-7/25/(exp(1/5*kmrgd(1)+27/5)-1)+1/5*(7/25*kmrgd(1)+189/25)/(exp(1/5*kmrgd(1)+27/5)-1)^2*exp(1/5*kmrgd(1)+27/5),0,0,0,0];[-8/625*(1-kmrgd(3))/(1-exp(-1/5*kmrgd(1)-52/5))^2*exp(-1/5*kmrgd(1)-52/5)+2/25*(1-kmrgd(3))*(4/125*kmrgd(1)+208/125)/(1-exp(-1/5*kmrgd(1)-52/5))^3*exp(-1/5*kmrgd(1)-52/5)^2+1/25*(1-kmrgd(3))*(4/125*kmrgd(1)+208/125)/(1-exp(-1/5*kmrgd(1)-52/5))^2*exp(-1/5*kmrgd(1)-52/5)-1/3200*kmrgd(3)*exp(-57/40-1/40*kmrgd(1)),0,-4/125/(1-exp(-1/5*kmrgd(1)-52/5))+1/5*(4/125*kmrgd(1)+208/125)/(1-exp(-1/5*kmrgd(1)-52/5))^2*exp(-1/5*kmrgd(1)-52/5)+1/80*exp(-57/40-1/40*kmrgd(1)),0,0,0];[4/10125*(1-kmrgd(4))*exp(-25/9-1/18*kmrgd(1))-8/25*kmrgd(4)/(1+exp(-1/5*kmrgd(1)-27/5))^3*exp(-1/5*kmrgd(1)-27/5)^2+4/25*kmrgd(4)/(1+exp(-1/5*kmrgd(1)-27/5))^2*exp(-1/5*kmrgd(1)-27/5),0,0,8/1125*exp(-25/9-1/18*kmrgd(1))-4/5/(1+exp(-1/5*kmrgd(1)-27/5))^2*exp(-1/5*kmrgd(1)-27/5),0,0];[1/5000/(1+exp(-1/10*kmrgd(1)-7/2))^3*exp(-1/10*kmrgd(1)-7/2)^2*(33/10*exp(1/20*kmrgd(1)+7/4)+exp(-1/20*kmrgd(1)-7/4))-1/10000/(1+exp(-1/10*kmrgd(1)-7/2))^2*exp(-1/10*kmrgd(1)-7/2)*(33/10*exp(1/20*kmrgd(1)+7/4)+exp(-1/20*kmrgd(1)-7/4))+1/500/(1+exp(-1/10*kmrgd(1)-7/2))^2*exp(-1/10*kmrgd(1)-7/2)*(33/200*exp(1/20*kmrgd(1)+7/4)-1/20*exp(-1/20*kmrgd(1)-7/4))+1/100*(1/(1+exp(-1/10*kmrgd(1)-7/2))-kmrgd(5))*(33/4000*exp(1/20*kmrgd(1)+7/4)+1/400*exp(-1/20*kmrgd(1)-7/4)),0,0,0,-33/20000*exp(1/20*kmrgd(1)+7/4)+1/2000*exp(-1/20*kmrgd(1)-7/4),0];[-2/3125*gca/(1+exp(-2/5*kmrgd(1)-10))^3*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)^2-1/625*gca/(1+exp(-2/5*kmrgd(1)-10))^2*exp(-2/5*kmrgd(1)-10)+1/3125*gca/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10),0,0,0,0,0]];
hess2=[[-3*gna*kmrgd(2)^2*kmrgd(4)/-C,-6*gna*kmrgd(2)*kmrgd(4)*(kmrgd(1)-ena)/-C,0,-3*gna*kmrgd(2)^2*(kmrgd(1)-ena)/-C,0,0];[-8/25/(1-exp(-27/2-1/4*kmrgd(1)))+1/4*(432/25+8/25*kmrgd(1))/(1-exp(-27/2-1/4*kmrgd(1)))^2*exp(-27/2-1/4*kmrgd(1))-7/25/(exp(1/5*kmrgd(1)+27/5)-1)+1/5*(7/25*kmrgd(1)+189/25)/(exp(1/5*kmrgd(1)+27/5)-1)^2*exp(1/5*kmrgd(1)+27/5),0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hess3=[[-4*gk*kmrgd(3)^3/-C,0,-12*gk*kmrgd(3)^2*(kmrgd(1)-ek)/-C,0,0,0];[0,0,0,0,0,0];[-4/125/(1-exp(-1/5*kmrgd(1)-52/5))+1/5*(4/125*kmrgd(1)+208/125)/(1-exp(-1/5*kmrgd(1)-52/5))^2*exp(-1/5*kmrgd(1)-52/5)+1/80*exp(-57/40-1/40*kmrgd(1)),0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hess4=[[-gna*kmrgd(2)^3/-C,-3*gna*kmrgd(2)^2*(kmrgd(1)-ena)/-C,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[8/1125*exp(-25/9-1/18*kmrgd(1))-4/5/(1+exp(-1/5*kmrgd(1)-27/5))^2*exp(-1/5*kmrgd(1)-27/5),0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hess5=[[-gm/-C,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[-33/20000*exp(1/20*kmrgd(1)+7/4)+1/2000*exp(-1/20*kmrgd(1)-7/4),0,0,0,0,0];[0,0,0,0,0,0]];
hess6=[[(-gahp/(1+kmrgd(6))+gahp*kmrgd(6)/(1+kmrgd(6))^2)/-C,0,0,0,0,-(-2*gahp/(1+kmrgd(6))^2+2*gahp*kmrgd(6)/(1+kmrgd(6))^3)*(kmrgd(1)-ek)/-C];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
hess(:,:,5) =hess5;
hess(:,:,6) =hess6;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
hessp1=[[-(-gna*kmrgd(2)^3*kmrgd(4)-gk*kmrgd(3)^4-gm*kmrgd(5)-gahp*kmrgd(6)/(1+kmrgd(6))-gl-2/5*gca/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)-gca/(1+exp(-2/5*kmrgd(1)-10)))/-C^2,3*gna*kmrgd(2)^2*kmrgd(4)*(kmrgd(1)-ena)/-C^2,4*gk*kmrgd(3)^3*(kmrgd(1)-ek)/-C^2,gna*kmrgd(2)^3*(kmrgd(1)-ena)/-C^2,gm*(kmrgd(1)-ek)/-C^2,(gahp/(1+kmrgd(6))-gahp*kmrgd(6)/(1+kmrgd(6))^2)*(kmrgd(1)-ek)/-C^2];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp2=[[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp3=[[-kmrgd(2)^3*kmrgd(4)/-C,-3*kmrgd(2)^2*kmrgd(4)*(kmrgd(1)-ena)/-C,0,-kmrgd(2)^3*(kmrgd(1)-ena)/-C,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp4=[[0,3*gna*kmrgd(2)^2*kmrgd(4)/-C,0,gna*kmrgd(2)^3/-C,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp5=[[-kmrgd(3)^4/-C,0,-4*kmrgd(3)^3*(kmrgd(1)-ek)/-C,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp6=[[0,0,4*gk*kmrgd(3)^3/-C,0,gm/-C,(gahp/(1+kmrgd(6))-gahp*kmrgd(6)/(1+kmrgd(6))^2)/-C];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp7=[[-1/-C,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp8=[[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp9=[[-kmrgd(5)/-C,0,0,0,-(kmrgd(1)-ek)/-C,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp10=[[-kmrgd(6)/(1+kmrgd(6))/-C,0,0,0,0,-(1/(1+kmrgd(6))-kmrgd(6)/(1+kmrgd(6))^2)*(kmrgd(1)-ek)/-C];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0]];
hessp11=[[(-2/5/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)-1/(1+exp(-2/5*kmrgd(1)-10)))/-C,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[-1/1250/(1+exp(-2/5*kmrgd(1)-10))^2*(kmrgd(1)-eca)*exp(-2/5*kmrgd(1)-10)-1/500/(1+exp(-2/5*kmrgd(1)-10)),0,0,0,0,0]];
hessp12=[[2/5*gca/(1+exp(-2/5*kmrgd(1)-10))^2*exp(-2/5*kmrgd(1)-10)/-C,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[0,0,0,0,0,0];[1/1250*gca/(1+exp(-2/5*kmrgd(1)-10))^2*exp(-2/5*kmrgd(1)-10),0,0,0,0,0]];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
hessp(:,:,10) =hessp10;
hessp(:,:,11) =hessp11;
hessp(:,:,12) =hessp12;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,C,I,gna,ena,gk,ek,gl,el,gm,gahp,gca,eca)
