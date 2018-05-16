% This function fills the matrix Q as defined in the text
% (see Appendix VII-B and C, equation (19) and (20))
% and returns the product Sp^4*A*Q

function M=matrixNparam(t,z,N)

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);
    z3(i)=z3(i-1)+cos(z3(2*N+i-1));
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1));
end

M3=matrix3Nparameters(t,z3,N);

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

end

        
