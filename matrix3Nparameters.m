% This function fills the matrix defined as A in the text 
% (see Appendix VII-C, equation (20) and following)
% and returns Sp^4 * A

function M=matrix3Nparameters(t,z,N)
global gamma Sp

x=z(1:N);
y=z(N+1:2*N);
th=z(2*N+1:3*N);
F=zeros(2,3*N);
T=zeros(N,3*N);

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

M=[F;T];
end


        
