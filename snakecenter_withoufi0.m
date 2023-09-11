function [xDot]=snakecenter_withoufi0(t,x,param)
%robot parameters
ct=param.ct;
cn=param.cn;
m=param.m;
l=param.l;
N=param.N;
J=m*l^2/3;
%control parameters
Kp=param.Kp;
Kd=param.Kd;
%state vector
fi=x(1:N);
p=x(N+1:N+2);
fiDot=x(N+3:2*N+2);
pDot=x(2*N+3:2*N+4);
%%%%%%%%%
H=ones(N,N);
H=triu(H);
theta=H*fi;
thetaDot=H*fiDot;
A=[diag(ones(N-1,1)) zeros(N-1,1)]+[zeros(N-1,1) diag(ones(N-1,1))];
D=[diag(ones(N-1,1)) zeros(N-1,1)]-[zeros(N-1,1) diag(ones(N-1,1))];
e=ones(N,1);
E=[e,zeros(N,1);zeros(N,1) e];
Stheta=diag(sin(theta));
Ctheta=diag(cos(theta));
V=A'*(D*D')^(-1)*A;
K=A'*(D*D')^(-1)*D;
W=m*l^2*Stheta*V*Ctheta-m*l^2*Ctheta*V*Stheta;
M=J*eye(N)+m*l^2*Stheta*V*Stheta+m*l^2*Ctheta*V*Ctheta;
u=zeros(N-1,1);
T1=[D'*(D*D')^(-1),e];
X=(T1)*[-l*A*cos(theta);p(1)];
Y=(T1)*[-l*A*sin(theta);p(2)];

% global fi0_last;
% global flag;
% 
% for i=1:size(flag,2)
%     %robot has not reached the point
%     if ~flag(i)
%         rx = param.Px(i)-p(1);
%         ry = param.Py(i)-p(2);
%         r  = sqrt(rx^2+ry^2);
%         if r>0.05
%             theta_bar = (sum(theta)/N);
%             theta_ref = -atan2(ry,rx);
%             fi0  = param.st1*(theta_bar+theta_ref);
%              if fi0 > deg2rad(40)
%                  fi0=deg2rad(40);
%              elseif fi0 < -deg2rad(40)
%                  fi0=-deg2rad(40);
%              end
%             fi0_last=fi0;
%         else
%             flag(i)=t;
%             fi0 = fi0_last;
%         end
%         break
%     else
%         fi0 = 0;
%     end
% end
for i = 1:N-1
%     fi_z = param.ampli_u*sin(param.omega_u*t+(i-1)*param.delta_u); 
%     fi_z = fi_z;%+fi0;
    fi_z = novel_wave(t, param.ampli_u, param.omega_u, i, param.delta_u);
    u(i) = Kp*(fi_z-fi(i))-Kd*fiDot(i); %control law
end

XDot=l*K'*Stheta*thetaDot+e*pDot(1);
YDot = -l*K'*Ctheta*thetaDot+e*pDot(2);
% equation 2.25
fR  = -[ct*(Ctheta.^2)+cn*(Stheta.^2) (ct-cn)*Stheta*Ctheta; ...
    (ct-cn)*Stheta*Ctheta ct*(Stheta.^2)+cn*(Ctheta.^2)]*[XDot;YDot];
% equation 2.39 in two parts - 
% 2.41a i 2.41b - actuated i unactuated joints
MFi = [H'*M*H, zeros(N,2); ...
      zeros(2,N), N*m*eye(2,2)];
% MFi_11 = MFi(1:N-1, 1:N-1);
% MFi_12 = MFi(N:N+2, 1:3);
MFi_21 = MFi(N:N+2, 1:N-1);
MFi_22 = MFi(N:N+2, N:N+2);
WFi = [H'*W*diag(H*fiDot)*H*fiDot; zeros(2,1)];
% WFi_1 = WFi(1:N-1);
WFi_2 = WFi(N:N+2);
% equation 2.40d
GFi = [-l*H'*Stheta*K, l*H'*Ctheta*K; ...
      -e', zeros(1,N); ...
      zeros(1,N), -e'];
% GFi_1 = GFi(1:N-1, 1:2*N);
GFi_2 = GFi(N:N+2, 1:2*N);
% equation 2.46a i 2.46b
AqFi = -MFi_22^(-1)*(WFi_2 +GFi_2*fR);
BqA  = -MFi_22^(-1)*MFi_21;
% State vector from equation 2.47
xDot=[fiDot; pDot; u; AqFi+BqA*u]; 
end