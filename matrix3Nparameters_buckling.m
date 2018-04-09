function M=matrix3Nparameters_buckling(t,z,N)
global gamma Sp LL

x=z(1:N);
y=z(N+1:2*N);
th=z(2*N+1:3*N);
F=zeros(2,3*N);
T=zeros(N,3*N);
C=zeros(2*(N-1),3*N);

for i=1:N
    u=cos(th(i));
    v=sin(th(i));
    F(1,i)=-(gamma*u^2+v^2);
    F(1,N+i)=-(gamma-1)*u*v;
    F(2,i)=-(gamma-1)*u*v;
    F(2,N+i)=-(u^2+gamma*v^2);
    F(1,2*N+i)=1/2*v;
    F(2,2*N+i)=-1/2*u;
end

F=Sp^3*F;

for i=1:N
    for j=i:N
        u=cos(th(j));
        v=sin(th(j));
        A=(x(j)-x(i));
        B=(y(j)-y(i));
        T(i,j)=1/2*v+...
            A*(-gamma+1)*v*u+...
            B*(gamma*u*u+v*v);
        T(i,N+j)=-1/2*u+...
            B*(gamma-1)*v*u-...
            A*(u*u+gamma*v*v);
        T(i,2*N+j)=-1/3-...
            1/2*A*u...
            -1/2*B*v;
    end
end

T=Sp^3*T;

for i=1:N-1
    C(i,i)=-1;
    C(N-1+i,N+i)=-1;
    C(i,2*N+i)=sin(th(i));
    C(N-1+i,2*N+i)=-cos(th(i));
    C(i,i+1)=1;
    C(N-1+i,N+i+1)=1;
end


M=[F;T;C];
M(1,:)=[-1,zeros(1,N-2),1,zeros(1,N),zeros(1,N-1),-sin(th(N))];
M(2,:)=[zeros(1,N),1,zeros(1,N-2),-1,zeros(1,N-1),-cos(th(N))];

end


        