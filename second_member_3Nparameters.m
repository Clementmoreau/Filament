function B=second_member_3Nparameters(t,z,N)
global Hx Hy
Mag=zeros(1,N+3);
th=[0;z(2*N+1:3*N)];

B1=zeros(3*N,1);
B2=zeros(3*N,1);
B3=zeros(3*N,1);
B1(N+2)=(th(N+1)-th(N));
B2(N+2)=-Mag(N)*sin(th(N+1));
B3(N+2)=Mag(N)*cos(th(N+1));
for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
    B2(3+i)=B2(3+i+1)-Mag(i+1)*sin(th(i+2));
    B3(3+i)=B3(3+i+1)+Mag(i+1)*cos(th(i+2));
end
B1(3)=cos(t);

BB=B1;%+Hx(t)*B2+Hy(t)*B3;
M=matrix3Nparameters(t,z,N);
B=M\BB;
end


    
    