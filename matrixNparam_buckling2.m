function M=matrixNparam_buckling2(t,z,N,Sp)

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);
    z3(i)=z3(i-1)+cos(z3(2*N+i-1));
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1));
end

M3=matrix3Nparameters_buckling2(t,z3,N,Sp);

C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;
    end
end
for i=2:N
    C1(i,i-1)=-sin(z3(2*N+i-1));
    C2(i,i-1)=cos(z3(2*N+i-1));
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-sin(z3(2*N+j));
        C2(i,j)=C2(i,j+1)+cos(z3(2*N+j));
    end
end
C=[ones(N,1),zeros(N,1),C1;zeros(N,1),ones(N,1),C2;zeros(N,2),C3];
M=M3*C;
U=zeros(N,1);
V=zeros(N,1);
U(N)=-sin(z3(3*N));
V(N)=+cos(z3(3*N));
for k=N-1:-1:1
    U(k)=U(k+1)-sin(z3(2*N+k));
    V(k)=V(k+1)+cos(z3(2*N+k));
end
%add 4 rows and 4 colums with unknown reactions
M1=[1,0,1,0;0,1,0,1;[zeros(N,2),U,V]];
%constraints rows
M2=[M1',zeros(4,4)];
% M2=[1,zeros(1,N+5);0,1,zeros(1,N+4);1,zeros(1,N+5);0,1,zeros(1,N+4)];
% M(3,N+2)=-sin(z3(3*N));
% M(4,N+2)=cos(z3(3*N));
% for i=N-1:-1:1
%     M2(3,i+2)=M2(3,i+3)-sin(z3(2*N+i));
%     M2(4,i+2)=M2(4,i+3)+cos(z3(2*N+i));
% end
M=[M,M1;M2];
end

        