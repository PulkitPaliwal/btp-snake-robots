clc;
clear;
close all;
param.N  = 10;
param.m  = 1;
param.l  = 0.07;
param.ct = 1;
param.cn = 10;
param.Kp  =7; %%%***
param.Kd  = 2; %%%***
param.st1 = 1; % path following parameter - for 0 - without following
param.st2 = 2; % path following parameter - can't be equal to zero (division)
param.st3 = 0.5; %%%***
t = 0:0.1:300;
param.delta_u = deg2rad(30);
param.ampli_u = deg2rad(20);
param.omega_u = deg2rad(50);

Px=[1; 3; 4; 6; 8];
Py=[1; 1; 0; -1; -1];
param.Px=Px;
param.Py=Py;
global flag;
flag=zeros(1,size(param.Px,1));
fi_r=zeros(1,param.N-1);

for i = 1:param.N-1
     fi_z = param.ampli_u*sin(param.omega_u*t+(i-1)*param.delta_u);
%      fi_z=2*i;
     fi_zad{i}=fi_z;
     fi_r(i)=fi_z(1);
end
%
theta    = zeros(param.N,1);
thetaDot = zeros(param.N,1);
fi       = fi_r;
fiDot    = zeros(param.N-1,1);
xCm     = 0;
yCm     = 0;
p       = [xCm; yCm];
xCmDot  = 0;
yCmDot  = 0;
pDot    = [xCmDot; yCmDot];
qA0    = fi';
qU0    = [theta(param.N); p];
qADot0 = fiDot;
qUDot0 = [thetaDot(param.N); pDot];
x0     = [qA0; qU0; qADot0; qUDot0];

[T1,X1]=ode45(@(t,x) snakecenter_withoufi0(t,x,param), t, x0); %undisturbed model
f1=flag;
info{1}=('Center tracking  SJC');
T={T1};
X={X1};
l=cell(param.N-1,1);
for j=1:size(X,2)
    Xa=X{j};
    Ta=T{j};
    figure
    hold on;
    grid on;
    
        
    for i = 1:param.N-1
        subplot(3,3,i)
        plot(Ta,Xa(:,i),'linewidth',3);
        hold on
%     plot(Ta,fi_zad(:,i));
   plot(Ta,fi_zad{i});
        ylabel(strcat('\phi_{', num2str(i), '}'));
        xlabel('t')
       %ylim([0,max(Xa(:,i))])
        xlim([0, max(Ta)])
        set(gca,'FontSize',12)
        grid on
        hold on;
         %plot(T1,fi_z(:,i));
        hold on
        plot(Ta,ones(size((Ta)))*deg2rad(40),'r:','linewidth',2);
        plot(Ta,ones(size((Ta)))*deg2rad(-40),'r:','linewidth',2);
        hold off
        
              if i==2
            title(info{j},'FontSize',12)
        end
    end
    
    
end
figure()
Xa=X{1};
plot(Xa(:,param.N+1),Xa(:,param.N+2),'linewidth',3);
hold on;
%  
 % plot(param.Px,param.Py,'o','MarkerSize',20,'Linewidth',4)

set(gca,'XLim',[min(param.Px)-1,max(param.Px)+1])
set(gca,'YLim',[min(-param.Py)-1,max(-param.Py)+1])
set(gca,'FontSize',15)
grid on
xlabel('x');
ylabel('y');
title('SJC for center tracking','FontSize',20)
A=[diag(ones(param.N-1,1)) zeros(param.N-1,1)]+[ zeros(param.N-1,1),diag(ones(param.N-1,1))];
D=[diag(ones(param.N-1,1)) zeros(param.N-1,1)]-[ zeros(param.N-1,1),diag(ones(param.N-1,1))];
e=ones(param.N,1);
T1=[D;e'./param.N];
Xpa=zeros(size(Xa,1),param.N);
Ypa=zeros(size(Xa,1),param.N);
for i=1:size(Xa,1)  % loop for time
    thetaa=zeros(1,param.N);
    
    
    for j=param.N:-1:1
        if j==param.N
            thetaa(param.N)=Xa(i,param.N);
            
        else
            thetaa(j)=Xa(i,j)+thetaa(j+1);
            %thetab(j)=Xb(i,j)+thetab(j+1)
           
        end
    end
    Tha(:,i)=thetaa';
    %Thb(:,i)=thetab';
    theetaabbb=Tha';
    pXa=Xa(i,param.N+1);
    pYa=Xa(i,param.N+2);
   
    %pXb=Xb(i,param.N+1);
    %pYb=Xb(i,param.N+2);
    l=0.07;
    Xpa(i,:)=T1^(-1)*[-A*cos(thetaa)'; pXa];
    Ypa(i,:)=T1^(-1)*[-A*sin(thetaa)'; pYa];
end
figure()
for i=1:length(Xpa)
clf;
j1=Xpa(i,:)+l*cos(Tha(:,i)');
j3=Ypa(i,:)+l*sin(Tha(:,i)');

plot(j1',j3','-o')
hold on

pp=plot(Xa(i,11),Xa(i,12),'o');
pp.MarkerFaceColor = [1 0.5 0];
ylim([-20,20])
xlim([-20,20])
grid on
hold on
movieframe4(i)=getframe;
end
% plot(j1',j3','-o')
% hold on
% 
% pp=plot(Xa(i,11),Xa(i,12),'o');
% pp.MarkerFaceColor = [1 0.5 0];
% ylim([-30,30])
% xlim([-30,30])
% 
% grid on
% hold on
% mywriter4=VideoWriter('snake5','MPEG-4');
% mywriter4.FrameRate=20;
% open(mywriter4);
% writeVideo(mywriter4,movieframe4);
% close(mywriter4);
% plot(j1(3001,:)',j3(3001,:)','-o')
% hold on
% pp=plot(Xa(i,11),-Xa(i,12),'o');
%  pp.MarkerFaceColor = [1 0.5 0];
%  ylim([-30,30]);
% xlim([-30,30]) ;
% grid on 
% figure 
% plot(j1(1,:)',j3(1,:)','-o')
% ylim([-30,30])
% xlim([-30,30])
% grid on
% hold on
% plot(X1(1,11),X1(1,12),'o')
% ylim([-30,30])
% xlim([-30,30])
% grid onx
% hold on
% plot(j1(100,:)',j3(100,:)','-o')
% ylim([-30,30])
% xlim([-30,30])
% grid on
%  hold on
% plot(X1(100,11),X1(100,12),'o')
% ylim([-30,30])
% xlim([-30,30])
% grid on
% hold on
% plot(j1(30001,:)',j3(30001,:)','-o')
% hold on
% plot(X1(30001,11),X1(30001,12),'o')
% grid on
% ylim([-30,30])
% xlim([-30,30])
