function S=second_member_bundle(t,z,N)

global gamma Sp k a th_max

% k_base=1/10;

L=sqrt((z(3*N+1:4*N)-z(1:N)).^2+(z(4*N+1:5*N)-z(N+1:2*N)).^2);
% L_base1=sqrt((z(3*N+2)-z(1))^2+(z(4*N+2)-z(N+1))^2);
% L_base2=sqrt((z(3*N+1)-z(2))^2+(z(4*N+1)-z(N+2))^2);

Fx=-(1-a./L).*(z(3*N+1:4*N)-z(1:N));
Fy=-(1-a./L).*(z(4*N+1:5*N)-z(N+1:2*N));

% Fx_base1=-k_base*(1-sqrt(1+a^2)/L_base1)*(z(3*N+2)-z(1));
% Fy_base1=-k_base*(1-sqrt(1+a^2)/L_base1)*(z(4*N+2)-z(N+1));
% 
% Fx_base2=-k_base*(1-sqrt(1+a^2)/L_base2)*(z(3*N+1)-z(2));
% Fy_base2=-k_base*(1-sqrt(1+a^2)/L_base2)*(z(4*N+1)-z(N+2));

XA=zeros(N-1);
YA=zeros(N-1);
XB=zeros(N-1);
YB=zeros(N-1);

for i=1:N
    for j=i:N
        XA(i,j)=z(j)-z(i);
        YA(i,j)=z(N+j)-z(N+i);
        XB(i,j)=z(3*N+j)-z(3*N+i);
        YB(i,j)=z(4*N+j)-z(4*N+i);
    end
end

TA=k*sum(XA*Fy-YA*Fx,2);
TB=-k*sum(XB*Fy-YB*Fx,2);

% T_base1=-Fx_base1*(z(4*N+1)-z(N+2))+Fy_base1*(z(3*N+1)-z(2));
% T_base2=-Fx_base2*(z(4*N+2)-z(N+1))+Fy_base2*(z(3*N+2)-z(1));

%TB(1)=TB(1)+T_base1;
%TA(1)=TA(1)+T_base2;


P1=[second_member_bundle1(t,z(1:3*N),N);second_member_bundle1(t,z(3*N+1:end),N)];

% angle oscillation
th=th_max*sin(t);
dth=th_max*cos(t);


%forced sliding
% d=1/4;
% xd=d*cos(t)/k;

%TA(1)=th_max*cos(t);

% free ends
% P2=P1+k*[sum(Fx);sum(Fy);TA;zeros(2*N-2,1);-sum(Fx);-sum(Fy);TB;zeros(2*N-2,1)];
% angle oscillation
 P2=P1+[0;0;TA;zeros(2*N-2,1);0;0;TB;zeros(2*N-2,1)];
% P2=P1+[a/2*dth*cos(th);a/2*dth*sin(th);TA;zeros(2*N-2,1);-a/2*dth*cos(th);-a/2*dth*sin(th);TB;zeros(2*N-2,1)];
% P2=P1+[0;0;TA;zeros(2*N-2,1);...
%     -k*sum(Fx)-Fx_base1-Fx_base2;-k*sum(Fy)-Fy_base1-Fy_base2;TB;zeros(2*N-2,1)];

% P2(3)=P2(3)+50*sin(t*pi)/N;
% P2(3*N+3)=P2(3*N+3)+50*sin(t*pi)/N;

% moving end
 P2(3)=dth;
% P2(3*N+3)=dth;

% sliding
% P2 = P1+k*[xd;0;TA;zeros(2*N-2,1);-xd;0;TB;zeros(2*N-2,1)];

M=M_bundle(t,z,N);
S=M\P2;

end




