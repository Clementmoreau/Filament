function B=second_member_Nparam_buckling2(t,z,N)

global th0 f0

[X,Y,TH]=coordinates_swimmerN(z,N);
X=X*N;
Y=Y*N;

s0=sin(th0);

B1=zeros(N+2,1);
B2=zeros(N+2,1);

B1(4:N)=z(4:N);

for i=N:-1:1
    Fi=1-(sin(TH(i))-2*Y(i))/N/s0;
    for j=i:N
        u=cos(TH(j));
        A=(X(j)-X(i));
        Tj=(2*A+u)*(N-2*Y(j)/s0)/2/N-(3*A+2*u)*sin(TH(j))/s0/3/N;     
        B2(i+2)=B2(i+2)+Tj;  
    end
    B2(2)=B2(2)+Fi;
end
B2(2)=0;

BB=B1+f0*B2;
M=matrixNparam(t,z,N);
B=M\BB;

end
